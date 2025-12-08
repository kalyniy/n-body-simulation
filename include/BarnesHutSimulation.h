// BarnesHutSimulation.h
#pragma once
#include "SimulationAlgorithm.h"
#include "Particle.hpp"
#include "Octree.hpp"
#include <vector>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif

class BarnesHutSimulation : public SimulationAlgorithm {
public:
    explicit BarnesHutSimulation(float theta = 0.5f,
                                 int   leafBucketSize = 8,
                                 int   maxDepth = 32);

    void setTheta(float t) { theta_ = t; }

    void computeStep(std::vector<particle_t>& particles,
                     const SimParams& params) override;

private:
    void computeAccelerations_(std::vector<particle_t>& particles,
                               const SimParams& params);
    void integrate_(std::vector<particle_t>& particles,
                    const SimParams& params);
    
#ifdef USE_MPI
    void computeAccelerationsMPI_(std::vector<particle_t>& particles,
                                  const SimParams& params);
    void integrateMPI_(std::vector<particle_t>& particles,
                       const SimParams& params);

    // MPI-friendly broadcasted tree
    std::vector<Octree::MpiTreeNode> mpi_nodes_;
    std::vector<int>                 mpi_leafIndices_;

    void traverseMpiTreeAccumulate_(int nodeId,
                                    int i,
                                    float G,
                                    float eps2,
                                    float theta,
                                    vector3_t& acc,
                                    const std::vector<particle_t>& particles) const;
#endif

    float              theta_;
    Octree::BuildParams bp_;
};
