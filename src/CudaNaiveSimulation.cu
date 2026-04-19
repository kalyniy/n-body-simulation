#include "CudaNaiveSimulation.h"
#include "CudaMisc.h"
#include <cstdio>
#include <cuda_runtime.h>

// Tiled shared-memory force kernel
extern __shared__ float smem[];

__global__ void naiveAccelKernel(
    const float* __restrict__ px, const float* __restrict__ py, const float* __restrict__ pz,
    const float* __restrict__ m,
    float* __restrict__ ax, float* __restrict__ ay, float* __restrict__ az,
    int n, float G, float eps2)
{
    const int TILE = blockDim.x;
    float* s_px = smem;
    float* s_py = smem + TILE;
    float* s_pz = smem + 2*TILE;
    float* s_m  = smem + 3*TILE;

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    float my_px = 0, my_py = 0, my_pz = 0;
    if (i < n) { my_px = px[i]; my_py = py[i]; my_pz = pz[i]; }

    float a_x = 0, a_y = 0, a_z = 0;
    for (int t = 0; t < n; t += TILE) {
        int j = t + threadIdx.x;
        s_px[threadIdx.x] = j < n ? px[j] : 0;
        s_py[threadIdx.x] = j < n ? py[j] : 0;
        s_pz[threadIdx.x] = j < n ? pz[j] : 0;
        s_m[threadIdx.x]  = j < n ? m[j]  : 0;
        __syncthreads();

        if (i < n) {
            int end = min(TILE, n - t);
            for (int k = 0; k < end; ++k) {
                float rx = s_px[k] - my_px;
                float ry = s_py[k] - my_py;
                float rz = s_pz[k] - my_pz;
                
                float r2 = rx*rx + ry*ry + rz*rz + eps2;
                float inv = rsqrtf(r2);
                float inv3 = inv*inv*inv;
                float s = G * s_m[k] * inv3;
                a_x += rx*s;
                a_y += ry*s;
                a_z += rz*s;
            }
        }
        __syncthreads();
    }
    if (i < n) { ax[i] = a_x; ay[i] = a_y; az[i] = a_z; }
}

__global__ void integrateKernel(
    float* px, float* py, float* pz,
    float* vx, float* vy, float* vz,
    const float* ax, const float* ay, const float* az,
    int n, float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    vx[i] += ax[i]*dt;
    vy[i] += ay[i]*dt;
    vz[i] += az[i]*dt;

    px[i] += vx[i]*dt;
    py[i] += vy[i]*dt;
    pz[i] += vz[i]*dt;
}

// ---- Host implementation ----

CudaNaiveSimulation::CudaNaiveSimulation(int bs) : block_size_(bs) {
    cudaDeviceProp prop;
    cudaCheck(cudaGetDeviceProperties(&prop, 0), "cudaGetDeviceProperties failed in ctor()");
    std::printf("[CUDA] Using device 0: %s (%d SMs, %.0f MB)\n",
                prop.name, prop.multiProcessorCount, prop.totalGlobalMem/(1024.0f * 1024.0f));
}

CudaNaiveSimulation::~CudaNaiveSimulation() 
{ 
    freeDevice_();
}

void CudaNaiveSimulation::ensureDevice_(size_t n) {
    if (n <= d_cap_) return;
    freeDevice_();
    size_t bytes = n * sizeof(float);

    cudaCheck(cudaMalloc(&d_px_, bytes), "cuda malloc d_px_");
    cudaCheck(cudaMalloc(&d_py_, bytes), "cuda malloc d_py_");
    cudaCheck(cudaMalloc(&d_pz_, bytes), "cuda malloc d_pz_");

    cudaCheck(cudaMalloc(&d_vx_, bytes), "cuda malloc d_vx_");
    cudaCheck(cudaMalloc(&d_vy_, bytes), "cuda malloc d_vy_");
    cudaCheck(cudaMalloc(&d_vz_, bytes), "cuda malloc d_vz_");

    cudaCheck(cudaMalloc(&d_ax_, bytes), "cuda malloc d_ax_");
    cudaCheck(cudaMalloc(&d_ay_, bytes), "cuda malloc d_ay_");
    cudaCheck(cudaMalloc(&d_az_, bytes), "cuda malloc d_az_");

    cudaCheck(cudaMalloc(&d_mass_, bytes), "cuda malloc d_mass_");

    d_cap_ = n;
    device_valid_ = false;
}

void CudaNaiveSimulation::uploadToDevice_(const ParticleSoA& s) {
    size_t b = s.size() * sizeof(float);
    cudaCheck(cudaMemcpy(d_px_, s.pos_x.data(), b, cudaMemcpyHostToDevice), "cuda memcpy d_px_ (host to device)");
    cudaCheck(cudaMemcpy(d_py_, s.pos_y.data(), b, cudaMemcpyHostToDevice), "cuda memcpy d_py_ (host to device)");
    cudaCheck(cudaMemcpy(d_pz_, s.pos_z.data(), b, cudaMemcpyHostToDevice), "cuda memcpy d_pz_ (host to device)");
    cudaCheck(cudaMemcpy(d_vx_, s.vel_x.data(), b, cudaMemcpyHostToDevice), "cuda memcpy d_vx_ (host to device)");
    cudaCheck(cudaMemcpy(d_vy_, s.vel_y.data(), b, cudaMemcpyHostToDevice), "cuda memcpy d_vy_ (host to device)");
    cudaCheck(cudaMemcpy(d_vz_, s.vel_z.data(), b, cudaMemcpyHostToDevice), "cuda memcpy d_vz_ (host to device)");
    cudaCheck(cudaMemcpy(d_mass_, s.mass.data(), b, cudaMemcpyHostToDevice), "cuda memcpy d_mass_ (host to device)");
    device_valid_ = true;
    device_dirty_ = false;
}

void CudaNaiveSimulation::downloadFromDevice_(ParticleSoA& s) {
    size_t b = s.size() * sizeof(float);
    cudaCheck(cudaMemcpy(s.pos_x.data(), d_px_, b, cudaMemcpyDeviceToHost), "cuda memcpy d_px_ (device to host)");
    cudaCheck(cudaMemcpy(s.pos_y.data(), d_py_, b, cudaMemcpyDeviceToHost), "cuda memcpy d_py_ (device to host)");
    cudaCheck(cudaMemcpy(s.pos_z.data(), d_pz_, b, cudaMemcpyDeviceToHost), "cuda memcpy d_pz_ (device to host)");
    cudaCheck(cudaMemcpy(s.vel_x.data(), d_vx_, b, cudaMemcpyDeviceToHost), "cuda memcpy d_vx_ (device to host)");
    cudaCheck(cudaMemcpy(s.vel_y.data(), d_vy_, b, cudaMemcpyDeviceToHost), "cuda memcpy d_vy_ (device to host)");
    cudaCheck(cudaMemcpy(s.vel_z.data(), d_vz_, b, cudaMemcpyDeviceToHost), "cuda memcpy d_vz_ (device to host)");
    device_dirty_ = false;
}

void CudaNaiveSimulation::freeDevice_() {
    auto f = [](float*& p) { 
        if (p) {
            cudaFree(p);
            p = nullptr;
        } 
    };

    f(d_px_); f(d_py_); f(d_pz_);
    f(d_vx_); f(d_vy_); f(d_vz_);
    f(d_ax_); f(d_ay_); f(d_az_);
    f(d_mass_);

    d_cap_ = 0; device_valid_ = false;
}

void CudaNaiveSimulation::computeStep(ParticleSoA& particles, const SimParams& params)
{
    const int n = static_cast<int>(particles.size());
    if (n == 0) return;

    ensureDevice_(n);

    // Upload only on first step (or if host data changed)
    if (!device_valid_) {
        uploadToDevice_(particles);
    }

    const int grid = (n + block_size_ - 1) / block_size_;
    const size_t smem = 4 * block_size_ * sizeof(float);

    naiveAccelKernel<<<grid, block_size_, smem>>>(
        d_px_, d_py_, d_pz_, d_mass_,
        d_ax_, d_ay_, d_az_, n, params.G, params.min_r2);
    cudaCheck(cudaGetLastError(), "launch naiveAccelKernel");
    cudaCheck(cudaDeviceSynchronize(), "sync naiveAccelKernel");

    integrateKernel<<<grid, block_size_>>>(
        d_px_, d_py_, d_pz_, d_vx_, d_vy_, d_vz_,
        d_ax_, d_ay_, d_az_, n, params.dt);
    cudaCheck(cudaGetLastError(), "launch integrateKernel");
    cudaCheck(cudaDeviceSynchronize(), "sync integrateKernel");

    device_dirty_ = true;

    // Download every step so host SoA stays current for checkpointing.
    // In a production system you'd only download on checkpoint steps.
    downloadFromDevice_(particles);
}
