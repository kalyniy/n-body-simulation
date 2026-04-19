// BarnesHutSimulation.h
#pragma once
#include "SimulationAlgorithm.h"
#include "ParticleSoA.h"
#include "Octree.hpp"
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

class BarnesHutSimulation : public SimulationAlgorithm {
public:
    explicit BarnesHutSimulation(float theta = 0.5f,
                                 int   leafBucketSize = 8,
                                 int   maxDepth = 32);
    void setTheta(float t) { theta_ = t; }
    void computeStep(ParticleSoA& particles,
                     const SimParams& params) override;

private:
    void computeAccelerations_(ParticleSoA& particles,
                               const SimParams& params);
    void integrate_(ParticleSoA& p,
                    const SimParams& params);
    
#ifdef USE_MPI
    void computeAccelerationsMPI_(ParticleSoA& p, const SimParams& params);
    void integrateMPI_(ParticleSoA& p, const SimParams& params);

    std::vector<Octree::FlatTreeNode>     mpi_nodes_;
    std::vector<Octree::FlatLeafParticle> mpi_leaves_;
    std::vector<float>                     globalPosMass_;   // 4 floats per particle: x,y,z,m

    void traverseMpiTree_(int nodeId,
                          int i,
                          const float* px,
                          const float* py,
                          const float* pz,
                          float G,
                          float eps2,
                          float theta,
                          float& ax,
                          float& ay,
                          float& az) const;
#endif

    float              theta_;
    Octree::BuildParams bp_;
};
