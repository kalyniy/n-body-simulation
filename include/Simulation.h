#pragma once
#include <cstddef>
#include <vector>
#include <string>
#include "Particle.hpp"

struct SimParams
{
    float G = 1.0f;
    float dt = 0.05f;
    float min_r2 = 1e-8f; // softening
};

class NBodySimulation
{
public:
    explicit NBodySimulation(SimParams params = {}) : params_(params) {}

    std::vector<particle_t> &particles() { return particles_; }
    const std::vector<particle_t> &particles() const { return particles_; }

    void clear();
    void reserve(std::size_t n);
    void setupSolarSystem(int worldX = 600, int worldY = 600, int worldZ = 600);
    void generateRandom(std::size_t n, int worldX, int worldY, int worldZ, float minMass = 1.f, float maxMass = 1000.f);

    // I/O / datasets
    void loadHACC(const std::string &dir);

    // One physics step (no rendering)
    void step();

    // Config
    void setG(float g) { params_.G = g; }
    void setDt(float d) { params_.dt = d; }
    const SimParams &params() const { return params_; }

private:
    void computeAccelerations_();
    void integrate_();

private:
    SimParams params_{};
    std::vector<particle_t> particles_;
};
