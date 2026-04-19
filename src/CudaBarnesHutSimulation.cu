#include "CudaBarnesHutSimulation.h"
#include "CudaMisc.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include <cuda_runtime.h>

#ifdef USE_MPI
#include <mpi.h>
#endif


// Warp size is 32 on all current NVIDIA hardware. The paper's force kernel
// relies on lockstep warp execution throughout.
#define WARP_SIZE 32

// Maximum traversal stack depth per warp (paper caps tree depth; 64 is 2x
// safety headroom). Sized in shared memory as one stack per warp.
#define BH_STACK 64

// Full mask for __all_sync / __syncwarp.
#define FULL_MASK 0xffffffffu

struct GpuTreeNode
{
    float center[3];
    float half;
    float com[3];
    float mass;
    int   child[8];
    int   leafOffset;
    int   leafCount;
    int   isLeaf;
};

struct GpuLeafParticle
{
    float pos[3];
    float mass;
};

static_assert(sizeof(GpuTreeNode) == sizeof(Octree::FlatTreeNode), "layout mismatch");
static_assert(sizeof(GpuLeafParticle) == sizeof(Octree::FlatLeafParticle), "layout mismatch");


// ---- Constant-memory kernel parameters (paper section 6.3.1) ----
// Parameters never change during a time step. Passing them via __constant__
// memory avoids re-transmission on every kernel invocation and frees up
// registers versus per-call arguments.
__constant__ float c_G;
__constant__ float c_epssq;     // eps^2  (paper: "epssq")
__constant__ float c_itolsq;    // 1 / theta^2 (paper: "itolsq")
__constant__ float c_dt;
__constant__ int   c_nNodes;
__constant__ int   c_n;


__global__ void bhAccelKernel(
    const float* __restrict__ px,
    const float* __restrict__ py,
    const float* __restrict__ pz,
    float* __restrict__ ax,
    float* __restrict__ ay,
    float* __restrict__ az,
    const GpuTreeNode*     __restrict__ nodes,
    const GpuLeafParticle* __restrict__ leaves)
{
    // Paper section 6.3.6: one stack per warp lives in shared memory, and
    // exactly one "current node" slot per warp is broadcast to all lanes of
    // that warp. This replaces the original per-thread local-memory stack and
    // cuts global traffic for node fetches by up to 32x.
    extern __shared__ int s_mem[];

    const int warpsPerBlock = blockDim.x / WARP_SIZE;
    const int warpIdInBlock = threadIdx.x / WARP_SIZE;
    const int laneId        = threadIdx.x & (WARP_SIZE - 1);

    // Layout inside s_mem:
    //   [0 .. warpsPerBlock*BH_STACK)            : per-warp traversal stacks
    //   [warpsPerBlock*BH_STACK .. +warpsPerBlock) : per-warp stack pointers
    //   [ .. ]                                   : per-warp shared node slot
    //       (sizeof(GpuTreeNode)/4 ints per warp)
    int* const s_stacks = s_mem;
    int* const s_depth  = s_mem + warpsPerBlock * BH_STACK;
    int* const s_nodes  = s_mem + warpsPerBlock * BH_STACK + warpsPerBlock;

    int* const myStack  = s_stacks + warpIdInBlock * BH_STACK;
    int* const myDepth  = s_depth  + warpIdInBlock;
    GpuTreeNode* const myNodeBuf = reinterpret_cast<GpuTreeNode*>(
        s_nodes + warpIdInBlock * (sizeof(GpuTreeNode) / sizeof(int)));

    // Paper section 6.3.1: persistent threads / grid-stride loop. Launch a
    // fixed number of blocks tuned to the hardware and let each thread walk
    // the body array in strides.
    const int gridStride = blockDim.x * gridDim.x;

    for (int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < c_n;
         i += gridStride)
    {
        const float pxi = px[i];
        const float pyi = py[i];
        const float pzi = pz[i];

        float a_x = 0.f;
        float a_y = 0.f;
        float a_z = 0.f;

        // Initialize this warp's stack with the root node. Lane 0 does the
        // write; __syncwarp() publishes it to every lane.
        if (laneId == 0)
        {
            myStack[0] = 0;
            *myDepth   = 1;
        }
        __syncwarp(FULL_MASK);

        while (*myDepth > 0)
        {
            // Paper section 6.3.6: "we allow only one thread per warp to read
            // the pertinent data, and then that thread makes the data
            // available to the other threads in the same warp by caching the
            // data in shared memory."
            int id = -1;
            if (laneId == 0)
            {
                *myDepth = *myDepth - 1;
                id = myStack[*myDepth];
                if (id >= 0 && id < c_nNodes)
                {
                    *myNodeBuf = nodes[id];
                }
            }
            // Broadcast the node id to every lane so they all take the same
            // branch below. __shfl_sync is how the other 31 lanes learn what
            // lane 0 popped.
            id = __shfl_sync(FULL_MASK, id, 0);

            // threadfence_block + syncwarp: ensure the shared-memory write
            // from lane 0 is visible to the rest of the warp before they read.
            __threadfence_block();
            __syncwarp(FULL_MASK);

            if (id < 0 || id >= c_nNodes) continue;

            const GpuTreeNode& nd = *myNodeBuf;
            if (nd.mass <= 0.f) continue;

            const float rx  = nd.com[0] - pxi;
            const float ry  = nd.com[1] - pyi;
            const float rz  = nd.com[2] - pzi;
            const float dsq = rx * rx + ry * ry + rz * rz;

            if (nd.isLeaf)
            {
                // Direct-sum interaction with every particle in the bucket.
                if (nd.leafOffset >= 0)
                {
                    const int base  = nd.leafOffset;
                    const int count = nd.leafCount;
                    for (int k = 0; k < count; ++k)
                    {
                        const GpuLeafParticle lp = leaves[base + k];

                        const float dx = lp.pos[0] - pxi;
                        const float dy = lp.pos[1] - pyi;
                        const float dz = lp.pos[2] - pzi;

                        const float d2   = dx * dx + dy * dy + dz * dz + c_epssq;
                        const float inv  = rsqrtf(d2);
                        const float inv3 = inv * inv * inv;
                        const float s    = c_G * lp.mass * inv3;

                        a_x += dx * s;
                        a_y += dy * s;
                        a_z += dz * s;
                    }
                }
                continue;
            }

            // Multipole-acceptance criterion (paper's form, section 6.3.6).
            // Classic BH test is (s_width / d) < theta, which needs a sqrt
            // and a division. Squaring both sides and inverting theta gives
            //     dsq * itolsq >= sw*sw
            // with itolsq = 1/(theta*theta). No sqrt, no division.
            const float sw = nd.half + nd.half;
            const bool  farEnough = (dsq * c_itolsq) >= (sw * sw);

            // Paper section 6.3.6: "Determining the union's border can be
            // done extremely quickly using the __all thread-voting function."
            // If every lane in the warp agrees the node is far enough, we
            // approximate for all of them. Otherwise we descend for all of
            // them. This eliminates divergence at the MAC test.
            const int allFar = __all_sync(FULL_MASK, farEnough ? 1 : 0);

            if (allFar)
            {
                const float d2s  = dsq + c_epssq;
                const float inv  = rsqrtf(d2s);
                const float inv3 = inv * inv * inv;
                const float s    = c_G * nd.mass * inv3;

                a_x += rx * s;
                a_y += ry * s;
                a_z += rz * s;
                continue;
            }

            // Push children onto the warp's shared stack. Lane 0 does the
            // writes so pushes are serialized and the stack pointer stays
            // consistent. Paper section 6.3.6: "uses only one thread per warp
            // to control the iteration stack."
            //
            // Children are packed to the front of child[] by the exporter
            // (see Octree::exportNodeRecursive_, paper section 6.3.4), so as
            // soon as we see -1 the rest are -1 too. Early-termination here
            // reduces both divergence and memory traffic.
            if (laneId == 0)
            {
                int d = *myDepth;
                #pragma unroll
                for (int c = 0; c < 8; ++c)
                {
                    const int ch = nd.child[c];
                    if (ch == -1) break;
                    if (d < BH_STACK)
                    {
                        myStack[d] = ch;
                        ++d;
                    }
                }
                *myDepth = d;
            }
            __syncwarp(FULL_MASK);
        }

        ax[i] = a_x;
        ay[i] = a_y;
        az[i] = a_z;
    }
}


__global__ void bhIntegrateKernel(
    float* px,
    float* py,
    float* pz,
    float* vx,
    float* vy,
    float* vz,
    const float* ax,
    const float* ay,
    const float* az)
{
    const int gridStride = blockDim.x * gridDim.x;
    const float dt = c_dt;

    for (int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < c_n;
         i += gridStride)
    {
        vx[i] += ax[i] * dt;
        vy[i] += ay[i] * dt;
        vz[i] += az[i] * dt;

        px[i] += vx[i] * dt;
        py[i] += vy[i] * dt;
        pz[i] += vz[i] * dt;
    }
}


// ---- Host ----

CudaBarnesHutSimulation::CudaBarnesHutSimulation(float theta, int lb, int md, int bs)
    : theta_(theta), block_size_(bs)
{
    bp_.bucket_size = lb;
    bp_.max_depth   = md;
    bp_.bounds_pad  = 1e-2f;

    // The force kernel's shared-memory layout assumes block_size_ is a
    // multiple of the warp size.
    if (block_size_ <= 0 || (block_size_ % WARP_SIZE) != 0)
    {
        cudaCheck(cudaErrorInvalidValue,
                  "block_size_ must be a positive multiple of warp size (32)");
    }

#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int deviceCount;
    cudaCheck(cudaGetDeviceCount(&deviceCount), "cudaGetDeviceCount");

    int deviceId = rank % deviceCount;
    cudaCheck(cudaSetDevice(deviceId), "cudaSetDevice");

    cudaDeviceProp props;
    cudaCheck(cudaGetDeviceProperties(&props, deviceId), "cudaGetDeviceProperties");

    std::printf("[CUDA BH+MPI] Rank %d -> dev %d: %s\n", rank, deviceId, props.name);
#else
    cudaDeviceProp props;
    cudaCheck(cudaGetDeviceProperties(&props, 0), "cudaGetDeviceProperties");

    std::printf("[CUDA BH] Using device 0: %s (%d SMs)\n", props.name, props.multiProcessorCount);
#endif
}

CudaBarnesHutSimulation::~CudaBarnesHutSimulation()
{
    freeDevice_();
}

void CudaBarnesHutSimulation::ensureParticleDevice_(size_t n)
{
    if (n <= d_pcap_) return;

    auto f = [](float*& p) {
        if (p) {
            cudaFree(p);
            p = nullptr;
        }
    };

    f(d_px_);
    f(d_py_);
    f(d_pz_);
    f(d_vx_);
    f(d_vy_);
    f(d_vz_);
    f(d_ax_);
    f(d_ay_);
    f(d_az_);
    f(d_mass_);

    size_t b = n * sizeof(float);

    cudaCheck(cudaMalloc(&d_px_, b), "cudaMalloc d_px_");
    cudaCheck(cudaMalloc(&d_py_, b), "cudaMalloc d_py_");
    cudaCheck(cudaMalloc(&d_pz_, b), "cudaMalloc d_pz_");

    cudaCheck(cudaMalloc(&d_vx_, b), "cudaMalloc d_vx_");
    cudaCheck(cudaMalloc(&d_vy_, b), "cudaMalloc d_vy_");
    cudaCheck(cudaMalloc(&d_vz_, b), "cudaMalloc d_vz_");

    cudaCheck(cudaMalloc(&d_ax_, b), "cudaMalloc d_ax_");
    cudaCheck(cudaMalloc(&d_ay_, b), "cudaMalloc d_ay_");
    cudaCheck(cudaMalloc(&d_az_, b), "cudaMalloc d_az_");

    cudaCheck(cudaMalloc(&d_mass_, b), "cudaMalloc d_mass_");

    d_pcap_       = n;
    device_valid_ = false;
}

void CudaBarnesHutSimulation::ensureTreeDevice_(size_t nN, size_t nL)
{
    if (nN > d_ncap_)
    {
        if (d_nodes_) cudaFree(d_nodes_);
        cudaCheck(cudaMalloc(&d_nodes_, nN * sizeof(GpuTreeNode)), "cudaMalloc d_nodes_");
        d_ncap_ = nN;
    }

    if (nL > d_lcap_)
    {
        if (d_leaves_) cudaFree(d_leaves_);
        cudaCheck(cudaMalloc(&d_leaves_, nL * sizeof(GpuLeafParticle)), "cudaMalloc d_leaves_");
        d_lcap_ = nL;
    }
}

void CudaBarnesHutSimulation::freeDevice_()
{
    auto f = [](float*& p) {
        if (p) {
            cudaFree(p);
            p = nullptr;
        }
    };

    f(d_px_);
    f(d_py_);
    f(d_pz_);
    f(d_vx_);
    f(d_vy_);
    f(d_vz_);
    f(d_ax_);
    f(d_ay_);
    f(d_az_);
    f(d_mass_);

    if (d_nodes_)
    {
        cudaFree(d_nodes_);
        d_nodes_ = nullptr;
    }
    if (d_leaves_)
    {
        cudaFree(d_leaves_);
        d_leaves_ = nullptr;
    }

    d_pcap_ = 0;
    d_ncap_ = 0;
    d_lcap_ = 0;
}

void CudaBarnesHutSimulation::uploadParticles_(const ParticleSoA& s)
{
    size_t b = s.size() * sizeof(float);

    cudaCheck(cudaMemcpy(d_px_,   s.pos_x.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail pos_x (host to device)");
    cudaCheck(cudaMemcpy(d_py_,   s.pos_y.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail pos_y (host to device)");
    cudaCheck(cudaMemcpy(d_pz_,   s.pos_z.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail pos_z (host to device)");

    cudaCheck(cudaMemcpy(d_vx_,   s.vel_x.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail vel_x (host to device)");
    cudaCheck(cudaMemcpy(d_vy_,   s.vel_y.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail vel_y (host to device)");
    cudaCheck(cudaMemcpy(d_vz_,   s.vel_z.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail vel_z (host to device)");

    cudaCheck(cudaMemcpy(d_mass_, s.mass.data(),  b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail mass (host to device)");

    device_valid_ = true;
}

void CudaBarnesHutSimulation::downloadParticles_(ParticleSoA& s)
{
    size_t b = s.size() * sizeof(float);

    cudaCheck(cudaMemcpy(s.pos_x.data(), d_px_, b, cudaMemcpyDeviceToHost), "cudaMemcpy pos_x (device to host)");
    cudaCheck(cudaMemcpy(s.pos_y.data(), d_py_, b, cudaMemcpyDeviceToHost), "cudaMemcpy pos_y (device to host)");
    cudaCheck(cudaMemcpy(s.pos_z.data(), d_pz_, b, cudaMemcpyDeviceToHost), "cudaMemcpy pos_z (device to host)");

    cudaCheck(cudaMemcpy(s.vel_x.data(), d_vx_, b, cudaMemcpyDeviceToHost), "cudaMemcpy vel_x (device to host)");
    cudaCheck(cudaMemcpy(s.vel_y.data(), d_vy_, b, cudaMemcpyDeviceToHost), "cudaMemcpy vel_y (device to host)");
    cudaCheck(cudaMemcpy(s.vel_z.data(), d_vz_, b, cudaMemcpyDeviceToHost), "cudaMemcpy vel_z (device to host)");
}

void CudaBarnesHutSimulation::launchKernels_(int n, float G, float eps2, float theta, float dt)
{
    const int nN = static_cast<int>(flat_nodes_.size());
    const int nL = static_cast<int>(flat_leaves_.size());

    ensureTreeDevice_(nN, nL);

    cudaCheck(cudaMemcpy(d_nodes_, flat_nodes_.data(),
                         nN * sizeof(GpuTreeNode),
                         cudaMemcpyHostToDevice),
              "launchKernels_ cudaMemcpy fail d_nodes_ (host to device)");

    if (nL > 0)
    {
        cudaCheck(cudaMemcpy(d_leaves_, flat_leaves_.data(),
                             nL * sizeof(GpuLeafParticle),
                             cudaMemcpyHostToDevice),
                  "launchKernels_ cudaMemcpy fail d_leaves_ (host to device)");
    }

    // Push per-step parameters into constant memory. Paper section 6.3.1.
    const float itolsq = 1.f / (theta * theta);

    cudaCheck(cudaMemcpyToSymbol(c_G,      &G,      sizeof(float)), "cudaMemcpyToSymbol c_G");
    cudaCheck(cudaMemcpyToSymbol(c_epssq,  &eps2,   sizeof(float)), "cudaMemcpyToSymbol c_epssq");
    cudaCheck(cudaMemcpyToSymbol(c_itolsq, &itolsq, sizeof(float)), "cudaMemcpyToSymbol c_itolsq");
    cudaCheck(cudaMemcpyToSymbol(c_dt,     &dt,     sizeof(float)), "cudaMemcpyToSymbol c_dt");
    cudaCheck(cudaMemcpyToSymbol(c_nNodes, &nN,     sizeof(int)),   "cudaMemcpyToSymbol c_nNodes");
    cudaCheck(cudaMemcpyToSymbol(c_n,      &n,      sizeof(int)),   "cudaMemcpyToSymbol c_n");

    // Persistent-threads launch config (paper section 6.3.1). Pick a grid
    // sized to the hardware rather than to n, then let each thread stride
    // through the particle array.
    int deviceId;
    cudaCheck(cudaGetDevice(&deviceId), "cudaGetDevice");
    cudaDeviceProp props;
    cudaCheck(cudaGetDeviceProperties(&props, deviceId), "cudaGetDeviceProperties (launch)");

    const int warpsPerBlock = block_size_ / WARP_SIZE;

    // Shared-memory footprint for the force kernel:
    //   per-warp stack        : BH_STACK * sizeof(int)
    //   per-warp stack pointer: sizeof(int)
    //   per-warp node cache   : sizeof(GpuTreeNode)
    const size_t bytesPerWarp =
        BH_STACK * sizeof(int) +
        sizeof(int) +
        sizeof(GpuTreeNode);
    const size_t shmemBytes = warpsPerBlock * bytesPerWarp;

    // Aim for roughly 2 resident blocks per SM for the force kernel; fall
    // back to 1x SM count for the integrator which has almost no state.
    const int forceBlocks     = props.multiProcessorCount * 2;
    const int integrateBlocks = props.multiProcessorCount;

    bhAccelKernel<<<forceBlocks, block_size_, shmemBytes>>>(
        d_px_, d_py_, d_pz_,
        d_ax_, d_ay_, d_az_,
        (const GpuTreeNode*)d_nodes_,
        (const GpuLeafParticle*)d_leaves_);

    bhIntegrateKernel<<<integrateBlocks, block_size_>>>(
        d_px_, d_py_, d_pz_,
        d_vx_, d_vy_, d_vz_,
        d_ax_, d_ay_, d_az_);

    cudaCheck(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
}

void CudaBarnesHutSimulation::computeStep(ParticleSoA& p, const SimParams& params)
{
#ifdef USE_MPI
    computeStepMPI_(p, params);
#else
    computeStepSerial_(p, params);
#endif
}

void CudaBarnesHutSimulation::computeStepSerial_(ParticleSoA& p, const SimParams& params)
{
    int n = static_cast<int>(p.size());
    if (!n) return;

    // Download positions from GPU if device is ahead of host (for tree build on CPU)
    ensureParticleDevice_(n);
    if (!device_valid_)
    {
        uploadParticles_(p);
    }
    else
    {
        downloadParticles_(p);  // need current positions on host for tree build
    }

    // Build tree on CPU from host SoA
    Octree tree;
    tree.build(p.pos_x.data(), p.pos_y.data(), p.pos_z.data(), p.mass.data(), n, bp_);
    tree.exportToFlatTree(flat_nodes_, flat_leaves_);

    // GPU traversal + integration
    launchKernels_(n, params.G, params.min_r2, theta_, params.dt);

    // Download for checkpointing
    downloadParticles_(p);
}

#ifdef USE_MPI
void CudaBarnesHutSimulation::computeStepMPI_(ParticleSoA& p, const SimParams& params)
{
    int rank;
    int nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int localN = static_cast<int>(p.size());

    std::vector<int> counts(nproc);
    MPI_Allgather(&localN, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int globalN = 0;
    std::vector<int> displs(nproc, 0);
    if (rank == 0)
    {
        for (int r = 0; r < nproc; ++r)
        {
            displs[r] = globalN;
            globalN += counts[r];
        }
    }

    // Ensure host has current data
    ensureParticleDevice_(localN);
    if (device_valid_)
    {
        downloadParticles_(p);
    }
    else
    {
        uploadParticles_(p);
    }

    std::vector<float> sendbuf(4 * localN);
    for (int i = 0; i < localN; ++i)
    {
        sendbuf[4 * i]     = p.pos_x[i];
        sendbuf[4 * i + 1] = p.pos_y[i];
        sendbuf[4 * i + 2] = p.pos_z[i];
        sendbuf[4 * i + 3] = p.mass[i];
    }

    std::vector<int> rc;
    std::vector<int> rd;
    if (rank == 0)
    {
        rc.resize(nproc);
        rd.resize(nproc);

        int offset = 0;
        for (int r = 0; r < nproc; ++r)
        {
            rc[r] = 4 * counts[r];
            rd[r] = offset;
            offset += rc[r];
        }
        globalPosMass_.resize(offset);
    }

    MPI_Gatherv(sendbuf.data(),
                4 * localN,
                MPI_FLOAT,
                rank == 0 ? globalPosMass_.data() : nullptr,
                rank == 0 ? rc.data()             : nullptr,
                rank == 0 ? rd.data()             : nullptr,
                MPI_FLOAT,
                0,
                MPI_COMM_WORLD);

    int nodeCount = 0;
    int leafCount = 0;
    if (rank == 0)
    {
        std::vector<float> gpx(globalN);
        std::vector<float> gpy(globalN);
        std::vector<float> gpz(globalN);
        std::vector<float> gm(globalN);

        for (int i = 0; i < globalN; ++i)
        {
            gpx[i] = globalPosMass_[4 * i];
            gpy[i] = globalPosMass_[4 * i + 1];
            gpz[i] = globalPosMass_[4 * i + 2];
            gm[i]  = globalPosMass_[4 * i + 3];
        }

        Octree tree;
        tree.build(gpx.data(), gpy.data(), gpz.data(), gm.data(), globalN, bp_);
        tree.exportToFlatTree(flat_nodes_, flat_leaves_);

        nodeCount = flat_nodes_.size();
        leafCount = flat_leaves_.size();
    }

    MPI_Bcast(&nodeCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&leafCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    flat_nodes_.resize(nodeCount);
    flat_leaves_.resize(leafCount);

    if (nodeCount > 0)
    {
        MPI_Bcast(flat_nodes_.data(),
                  nodeCount * sizeof(Octree::FlatTreeNode),
                  MPI_BYTE,
                  0,
                  MPI_COMM_WORLD);
    }
    if (leafCount > 0)
    {
        MPI_Bcast(flat_leaves_.data(),
                  leafCount * sizeof(Octree::FlatLeafParticle),
                  MPI_BYTE,
                  0,
                  MPI_COMM_WORLD);
    }

    if (!device_valid_) uploadParticles_(p);
    launchKernels_(localN, params.G, params.min_r2, theta_, params.dt);
    downloadParticles_(p);
}
#endif
