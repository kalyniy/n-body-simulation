#include "BarnesHutSimulation.h"
#include <cmath>

void BarnesHutSimulation::computeStep(std::vector<particle_t>& particles, const SimParams& params)
{
    computeAccelerations_(particles, params);
    integrate_(particles, params);
}

void BarnesHutSimulation::computeAccelerations_(std::vector<particle_t>& particles, const SimParams& params)
{
    // Reset accelerations
    for (auto& p : particles) p.acceleration = {0, 0, 0};

    if (particles.empty()) return;

    // Build octree from current state
    Octree tree;
    tree.build(particles, bp_);

    const float G    = params.G;
    const float eps2 = params.min_r2;

    // Accumulate BH forces
    for (int i = 0; i < static_cast<int>(particles.size()); ++i)
    {
        // Acceleration from the tree (monopole + per-leaf exact)
        vector3_t ai = tree.accelerationOn(i, G, eps2, theta_);
        particles[i].acceleration = ai;
    }
}

void BarnesHutSimulation::integrate_(std::vector<particle_t>& particles, const SimParams& params)
{
    const float dt = params.dt;
    for (auto& p : particles) {
        p.velocity += p.acceleration * (0.5f * dt);  // half-kick
        p.position += p.velocity * dt;                // drift
        // [recompute accelerations]
        p.velocity += p.acceleration * (0.5f * dt); 
    }
}