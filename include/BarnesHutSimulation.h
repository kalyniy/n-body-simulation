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
    void computeStep(std::vector<particle_t>& particles, const SimParams& params) override;

private:
    void computeAccelerations_(std::vector<particle_t>& particles, const SimParams& params);
    void integrate_(std::vector<particle_t>& particles, const SimParams& params);
    
#ifdef USE_MPI
    void computeAccelerationsMPI_(std::vector<particle_t>& particles, const SimParams& params);
    void integrateMPI_(std::vector<particle_t>& particles, const SimParams& params);
    
    // MPI-specific helpers
    struct TreeNode {
        vector3_t com;
        float mass;
        float box_half;
        vector3_t box_center;
        bool is_leaf;
    };
    
    void extractEssentialNodes_(const Octree& tree, std::vector<TreeNode>& nodes);
    void buildLocalTree_(const std::vector<particle_t>& particles, 
                        size_t start, size_t end,
                        Octree& tree);
    void exchangeTreeData_(std::vector<TreeNode>& local_nodes,
                          std::vector<std::vector<TreeNode>>& all_nodes);
#endif

    float theta_;
    Octree::BuildParams bp_;
};