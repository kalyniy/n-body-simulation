#pragma once
#include "Particle.hpp"
#include <vector>
struct SimParams
{
    float G = 1.0f;
    float dt = 0.05f;
    float min_r2 = 1e-8f; // softening
};

class SimulationAlgorithm {
public:
    virtual ~SimulationAlgorithm() = default;
    virtual void computeStep(std::vector<particle_t>& particles, const SimParams& params) = 0;
};