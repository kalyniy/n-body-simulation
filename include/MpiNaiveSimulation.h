#pragma once
#include "SimulationAlgorithm.h"
#include "ParticleSoA.h"

class MpiNaiveSimulation : public SimulationAlgorithm {
public:
    void computeStep(ParticleSoA& particles, const SimParams& params) override;

private:
    void computeAccelerations_(ParticleSoA& particles, const SimParams& params);
    void integrate_(ParticleSoA& particles, const SimParams& params);
};