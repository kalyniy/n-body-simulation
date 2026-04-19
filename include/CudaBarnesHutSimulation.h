#pragma once

#include "SimulationAlgorithm.h"
#include "ParticleSoA.h"
#include "Octree.hpp"

#include <vector>

class CudaBarnesHutSimulation : public SimulationAlgorithm
{
public:
    explicit CudaBarnesHutSimulation(float theta = 0.5f,
                                     int leafBucket = 8,
                                     int maxDepth = 32,
                                     int blockSize = 256);

    ~CudaBarnesHutSimulation() override;

    void setTheta(float t) { theta_ = t; }

    void computeStep(ParticleSoA& particles, const SimParams& params) override;

private:
    void computeStepSerial_(ParticleSoA& p, const SimParams& params);

#ifdef USE_MPI
    void computeStepMPI_(ParticleSoA& p, const SimParams& params);
    std::vector<float> globalPosMass_;
#endif

    void ensureParticleDevice_(size_t n);
    void ensureTreeDevice_(size_t nN, size_t nL);
    void freeDevice_();

    void uploadParticles_(const ParticleSoA& s);
    void downloadParticles_(ParticleSoA& s);

    void launchKernels_(int n, float G, float eps2, float theta, float dt);

    float theta_;
    Octree::BuildParams bp_;
    int block_size_;
    bool device_valid_ = false;

    std::vector<Octree::FlatTreeNode>     flat_nodes_;
    std::vector<Octree::FlatLeafParticle> flat_leaves_;

    float* d_px_ = nullptr;
    float* d_py_ = nullptr;
    float* d_pz_ = nullptr;

    float* d_vx_ = nullptr;
    float* d_vy_ = nullptr;
    float* d_vz_ = nullptr;

    float* d_ax_ = nullptr;
    float* d_ay_ = nullptr;
    float* d_az_ = nullptr;

    float* d_mass_ = nullptr;
    size_t d_pcap_ = 0;

    void*  d_nodes_  = nullptr;
    void*  d_leaves_ = nullptr;
    size_t d_ncap_   = 0;
    size_t d_lcap_   = 0;
};