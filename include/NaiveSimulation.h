#pragma once
#include "SimulationAlgorithm.h"
#include "Particle.hpp"
#include "Simulation.h"

class NaiveSimulation : public SimulationAlgorithm {
public:
    void computeStep(std::vector<particle_t>& particles, const SimParams& params) override;

private:
    void computeAccelerations_(std::vector<particle_t>& particles, const SimParams& params);
    void integrate_(std::vector<particle_t>& particles, const SimParams& params);
};