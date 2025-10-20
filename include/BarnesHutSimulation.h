#pragma once
#include "SimulationAlgorithm.h"
#include "Particle.hpp"
#include "Octree.hpp"
#include <vector>
#include <iostream>
class BarnesHutSimulation : public SimulationAlgorithm {
public:
    explicit BarnesHutSimulation(float theta = 0.5f,
                                 int   leafBucketSize = 8,
                                 int   maxDepth = 32)
        : theta_(theta)
    {
        std::cout << "BarnesHutSimulation running\n";

        bp_.bucket_size = leafBucketSize;
        bp_.max_depth   = maxDepth;
        bp_.bounds_pad  = 1e-2f;
    }

    void setTheta(float t) { theta_ = t; }

    void computeStep(std::vector<particle_t>& particles, const SimParams& params) override;

private:
    void computeAccelerations_(std::vector<particle_t>& particles, const SimParams& params);
    void integrate_(std::vector<particle_t>& particles, const SimParams& params);

    float theta_;
    Octree::BuildParams bp_;
};