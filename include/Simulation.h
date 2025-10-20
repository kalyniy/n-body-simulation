#pragma once
#include <cstddef>
#include <vector>
#include <string>
#include <memory>
#include <atomic>
#include <mutex>
#include "Particle.hpp"
#include "SimulationAlgorithm.h"

class NBodySimulation
{
public:
    explicit NBodySimulation(std::unique_ptr<SimulationAlgorithm> algorithm, SimParams params = {})
        : algorithm_(std::move(algorithm)), params_(params) {}

    std::vector<particle_t> &particles() { return particles_; }
    const std::vector<particle_t> &particles() const { return particles_; }

    void clear();
    void reserve(std::size_t n);
    void setupSolarSystem(int worldX, int worldY, int worldZ);
    void generateRandom(std::size_t n, int worldX, int worldY, int worldZ, float minMass = 1.f, float maxMass = 1000.f);

    // I/O / datasets
    void loadHACC(const std::string &dir);

    // One physics step (no rendering)
    void step();
    const std::vector<particle_t>& getRenderBuffer();

    // Config
    void setG(float g) { params_.G = g; }
    void setDt(float d) { params_.dt = d; }
    void setAlgorithm(std::unique_ptr<SimulationAlgorithm> algorithm) { algorithm_ = std::move(algorithm); }

    const SimParams &params() const { return params_; }

private:
    std::unique_ptr<SimulationAlgorithm> algorithm_;
    SimParams params_{};
    std::vector<particle_t> particles_;
    std::vector<particle_t> render_buffer_;  // double buffer for rendering
    std::atomic<bool> buffer_ready_{false};
    std::mutex buffer_mutex_;
};
