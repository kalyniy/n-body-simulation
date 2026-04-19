#pragma once
#include "ParticleSoA.h"

struct SimParams
{
    float G = 1.0f;
    float dt = 0.05f;
    float min_r2 = 1e-8f; // softening
};

enum class AlgorithmKind {
    NaiveSeq,
    NaiveMpi,
    BarnesHutSeq,
    BarnesHutMpi,
    NaiveCuda,
    NaiveCudaMpi,
    BarnesHutCuda,
    BarnesHutCudaMpi
};

class SimulationAlgorithm {
public:
    virtual ~SimulationAlgorithm() = default;
    virtual void computeStep(ParticleSoA& particles, const SimParams& params) = 0;
};