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


#define BH_STACK 64


__global__ void bhAccelKernel(
    const float* __restrict__ px,
    const float* __restrict__ py,
    const float* __restrict__ pz,
    float* __restrict__ ax,
    float* __restrict__ ay,
    float* __restrict__ az,
    int n,
    const GpuTreeNode*     __restrict__ nodes,
    const GpuLeafParticle* __restrict__ leaves,
    int nNodes,
    float G,
    float eps2,
    float theta)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    float pxi = px[i];
    float pyi = py[i];
    float pzi = pz[i];

    float a_x = 0.f;
    float a_y = 0.f;
    float a_z = 0.f;

    int stack[BH_STACK];
    int sp = 0;
    stack[sp++] = 0;

    while (sp > 0)
    {
        int id = stack[--sp];
        if (id < 0 || id >= nNodes) continue;

        const auto& nd = nodes[id];
        if (nd.mass <= 0) continue;

        float rx = nd.com[0] - pxi;
        float ry = nd.com[1] - pyi;
        float rz = nd.com[2] - pzi;
        float r2 = rx * rx + ry * ry + rz * rz;

        if (nd.isLeaf)
        {
            if (nd.leafOffset >= 0)
            {
                for (int k = 0; k < nd.leafCount; ++k)
                {
                    const auto& lp = leaves[nd.leafOffset + k];

                    float dx = lp.pos[0] - pxi;
                    float dy = lp.pos[1] - pyi;
                    float dz = lp.pos[2] - pzi;

                    float d2 = dx * dx + dy * dy + dz * dz + eps2;
                    float inv = rsqrtf(d2);
                    float inv3 = inv * inv * inv;
                    float s = G * lp.mass * inv3;

                    a_x += dx * s;
                    a_y += dy * s;
                    a_z += dz * s;
                }
            }
            continue;
        }

        float sw = nd.half * 2.f;
        float d = sqrtf(r2) + 1e-12f;

        if ((sw / d) < theta)
        {
            float r2s = r2 + eps2;
            float inv = rsqrtf(r2s);
            float inv3 = inv * inv * inv;
            float s = G * nd.mass * inv3;

            a_x += rx * s;
            a_y += ry * s;
            a_z += rz * s;
            continue;
        }

        for (int c = 0; c < 8; ++c)
        {
            if (nd.child[c] != -1 && sp < BH_STACK)
            {
                stack[sp++] = nd.child[c];
            }
        }
    }

    ax[i] = a_x;
    ay[i] = a_y;
    az[i] = a_z;
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
    const float* az,
    int n,
    float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    vx[i] += ax[i] * dt;
    vy[i] += ay[i] * dt;
    vz[i] += az[i] * dt;

    px[i] += vx[i] * dt;
    py[i] += vy[i] * dt;
    pz[i] += vz[i] * dt;
}


// ---- Host ----

CudaBarnesHutSimulation::CudaBarnesHutSimulation(float theta, int lb, int md, int bs)
    : theta_(theta), block_size_(bs)
{
    bp_.bucket_size = lb;
    bp_.max_depth = md;
    bp_.bounds_pad = 1e-2f;

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

    d_pcap_ = n;
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

    cudaCheck(cudaMemcpy(d_px_, s.pos_x.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail pos_x (host to device)");
    cudaCheck(cudaMemcpy(d_py_, s.pos_y.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail pos_y (host to device)");
    cudaCheck(cudaMemcpy(d_pz_, s.pos_z.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail pos_z (host to device)");

    cudaCheck(cudaMemcpy(d_vx_, s.vel_x.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail vel_x (host to device)");
    cudaCheck(cudaMemcpy(d_vy_, s.vel_y.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail vel_y (host to device)");
    cudaCheck(cudaMemcpy(d_vz_, s.vel_z.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail vel_z (host to device)");

    cudaCheck(cudaMemcpy(d_mass_, s.mass.data(), b, cudaMemcpyHostToDevice), "uploadParticles_ cudaMemcpy fail mass (host to device)");

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
    int nN = flat_nodes_.size();
    int nL = flat_leaves_.size();

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

    int grid = (n + block_size_ - 1) / block_size_;

    bhAccelKernel<<<grid, block_size_>>>(
        d_px_, d_py_, d_pz_,
        d_ax_, d_ay_, d_az_,
        n,
        (GpuTreeNode*)d_nodes_,
        (GpuLeafParticle*)d_leaves_,
        nN,
        G, eps2, theta);

    bhIntegrateKernel<<<grid, block_size_>>>(
        d_px_, d_py_, d_pz_,
        d_vx_, d_vy_, d_vz_,
        d_ax_, d_ay_, d_az_,
        n,
        dt);

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
        sendbuf[4 * i] = p.pos_x[i];
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
            gm[i] = globalPosMass_[4 * i + 3];
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