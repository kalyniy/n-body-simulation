#pragma once
#include "SimulationAlgorithm.h"
#include "ParticleSoA.h"

class CudaNaiveSimulation : public SimulationAlgorithm
{
public:
    explicit CudaNaiveSimulation(int blockSize = 256);
    ~CudaNaiveSimulation() override;

    void computeStep(ParticleSoA& particles, const SimParams& params) override;

private:
    void ensureDevice_(size_t n);
    void uploadToDevice_(const ParticleSoA& soa);
    void downloadFromDevice_(ParticleSoA& soa);
    void freeDevice_();

    int block_size_;
    bool device_dirty_ = false;  // true = device has newer data than host
    bool device_valid_ = false;  // true = device buffers contain valid data

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
    size_t d_cap_ = 0;
};