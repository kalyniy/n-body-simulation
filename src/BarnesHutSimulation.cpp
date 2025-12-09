// BarnesHutSimulation.cpp
#include "BarnesHutSimulation.h"
#include <cmath>

BarnesHutSimulation::BarnesHutSimulation(float theta,
                                         int   leafBucketSize,
                                         int   maxDepth)
    : theta_(theta)
{
    std::printf("BarnesHutSimulation ctor()\n");
    bp_.bucket_size = leafBucketSize;
    bp_.max_depth   = maxDepth;
    bp_.bounds_pad  = 1e-2f;
}

void BarnesHutSimulation::computeStep(std::vector<particle_t>& particles,
                                      const SimParams& params)
{
#ifdef USE_MPI
    computeAccelerationsMPI_(particles, params);
    integrateMPI_(particles, params);
#else
    computeAccelerations_(particles, params);
    integrate_(particles, params);
#endif
}

// ======================= Serial (non-MPI) =======================

void BarnesHutSimulation::computeAccelerations_(std::vector<particle_t>& particles,
                                                const SimParams& params)
{
    for (auto& p : particles) {
        p.acceleration = {0, 0, 0};
    }
    if (particles.empty()) return;

    Octree tree;
    tree.build(particles, bp_);

    const float G    = params.G;
    const float eps2 = params.min_r2;

    for (int i = 0; i < static_cast<int>(particles.size()); ++i)
    {
        vector3_t ai = tree.accelerationOn(i, G, eps2, theta_);
        particles[i].acceleration = ai;
    }
}

void BarnesHutSimulation::integrate_(std::vector<particle_t>& particles,
                                     const SimParams& params)
{
    const float dt = params.dt;
    for (auto& p : particles) {
        p.velocity  += p.acceleration * dt;
        p.position  += p.velocity * dt;
    }
}

#ifdef USE_MPI

// ======================= MPI: Accelerations =======================
//
// Rank 0 builds the tree, flattens it, and broadcasts to all ranks.
// All ranks traverse the same tree but only for their subset of particles.

void BarnesHutSimulation::computeAccelerationsMPI_(std::vector<particle_t>& particles,
                                                   const SimParams& params)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int   localN = static_cast<int>(particles.size());
    const float G      = params.G;
    const float eps2   = params.min_r2;

    // Reset accelerations
    for (auto& p : particles) {
        p.acceleration = {0.0f, 0.0f, 0.0f};
    }
    if (localN == 0) return;

    // 1) Gather local counts so rank 0 knows the global N and layout
    std::vector<int> counts(size);
    MPI_Allgather(&localN, 1, MPI_INT,
                  counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    int globalN = 0;
    std::vector<int> displs(size, 0);
    if (rank == 0) {
        for (int r = 0; r < size; ++r) {
            displs[r] = globalN;
            globalN  += counts[r];
        }
    }

    // 2) Pack local positions + masses (4 floats per particle)
    std::vector<float> sendbuf(4 * localN);
    for (int i = 0; i < localN; ++i) {
        const auto& p = particles[i];
        const int idx = 4 * i;
        sendbuf[idx + 0] = p.position.x;
        sendbuf[idx + 1] = p.position.y;
        sendbuf[idx + 2] = p.position.z;
        sendbuf[idx + 3] = p.mass;
    }

    // 3) Gatherv packed data on rank 0 (like reference code's MPI_Gather on rr)
    std::vector<int> recvcountsF;
    std::vector<int> displsF;
    if (rank == 0) {
        recvcountsF.resize(size);
        displsF.resize(size);

        int offsetF = 0;
        for (int r = 0; r < size; ++r) {
            recvcountsF[r] = 4 * counts[r];
            displsF[r]     = offsetF;
            offsetF       += recvcountsF[r];
        }
        globalPosMass_.resize(offsetF);
    }

    MPI_Gatherv(sendbuf.data(), 4 * localN, MPI_FLOAT,
                rank == 0 ? globalPosMass_.data() : nullptr,
                rank == 0 ? recvcountsF.data()     : nullptr,
                rank == 0 ? displsF.data()         : nullptr,
                MPI_FLOAT,
                0, MPI_COMM_WORLD);

    // 4) Rank 0 builds tree from global positions/masses and serializes it
    int nodeCount = 0;
    int leafCount = 0;

    if (rank == 0) {
        globalParticles_.resize(globalN);

        for (int i = 0; i < globalN; ++i) {
            const int idx = 4 * i;
            auto& p = globalParticles_[i];

            p.position.x = globalPosMass_[idx + 0];
            p.position.y = globalPosMass_[idx + 1];
            p.position.z = globalPosMass_[idx + 2];
            p.mass       = globalPosMass_[idx + 3];

            // velocities/accels are not needed for tree build
            p.velocity    = {0.0f, 0.0f, 0.0f};
            p.acceleration = {0.0f, 0.0f, 0.0f};
        }

        Octree tree;
        tree.build(globalParticles_, bp_);
        tree.exportToMpiTree(mpi_nodes_, mpi_leafParticles_);

        nodeCount = static_cast<int>(mpi_nodes_.size());
        leafCount = static_cast<int>(mpi_leafParticles_.size());
    }

    // 5) Broadcast tree sizes
    MPI_Bcast(&nodeCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&leafCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    mpi_nodes_.resize(nodeCount);
    mpi_leafParticles_.resize(leafCount);

    // 6) Broadcast serialized tree
    if (nodeCount > 0) {
        MPI_Bcast(reinterpret_cast<unsigned char*>(mpi_nodes_.data()),
                  nodeCount * static_cast<int>(sizeof(Octree::MpiTreeNode)),
                  MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    if (leafCount > 0) {
        MPI_Bcast(reinterpret_cast<unsigned char*>(mpi_leafParticles_.data()),
                  leafCount * static_cast<int>(sizeof(Octree::MpiLeafParticle)),
                  MPI_BYTE, 0, MPI_COMM_WORLD);
    }

    // 7) Compute accelerations for this rank's *local* particles only
    for (int i = 0; i < localN; ++i) {
        vector3_t acc = {0.0f, 0.0f, 0.0f};

        if (nodeCount > 0) {
            traverseMpiTreeAccumulate_(0, i, G, eps2, theta_,
                                       acc, particles);
        }
        particles[i].acceleration = acc;
    }
}

// MPI tree traversal (mirrors Octree::traverseAccumulate_)

void BarnesHutSimulation::traverseMpiTreeAccumulate_(int nodeId,
                                                     int i,
                                                     float G,
                                                     float eps2,
                                                     float theta,
                                                     vector3_t& acc,
                                                     const std::vector<particle_t>& particles) const
{
    if (nodeId < 0 ||
        nodeId >= static_cast<int>(mpi_nodes_.size()))
    {
        return;
    }

    const auto& node = mpi_nodes_[nodeId];
    if (node.mass <= 0.0f) return;

    const vector3_t& pi = particles[i].position;

    vector3_t r = {
        node.com[0] - pi.x,
        node.com[1] - pi.y,
        node.com[2] - pi.z
    };
    float r2_true = r.x * r.x + r.y * r.y + r.z * r.z;

    if (node.isLeaf) {
        if (node.leafOffset >= 0 && node.leafCount > 0) {
            const int offset = node.leafOffset;
            const int count  = node.leafCount;

            for (int k = 0; k < count; ++k) {
                const auto& lp = mpi_leafParticles_[offset + k];

                vector3_t rij = {
                    lp.pos[0] - pi.x,
                    lp.pos[1] - pi.y,
                    lp.pos[2] - pi.z
                };

                float d2     = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z + eps2;
                float inv_r  = 1.0f / std::sqrt(d2);
                float inv_r3 = inv_r * inv_r * inv_r;
                float s      = G * lp.mass * inv_r3;

                acc.x += rij.x * s;
                acc.y += rij.y * s;
                acc.z += rij.z * s;
            }
        }
        return;
    }

    // Opening criterion: (s / d) < theta
    const float s = node.half * 2.0f;
    float d       = std::sqrt(r2_true) + 1e-12f;

    if ((s / d) < theta) {
        float r2_soft = r2_true + eps2;
        float inv_r   = 1.0f / std::sqrt(r2_soft);
        float inv_r3  = inv_r * inv_r * inv_r;
        float sgm     = G * node.mass * inv_r3;
        acc.x += r.x * sgm;
        acc.y += r.y * sgm;
        acc.z += r.z * sgm;
        return;
    }

    // Otherwise recurse into children
    for (int c = 0; c < 8; ++c) {
        int childId = node.child[c];
        if (childId != -1) {
            traverseMpiTreeAccumulate_(childId, i, G, eps2, theta,
                                       acc, particles);
        }
    }
}

// ======================= MPI: Integration =======================
//
// Each rank updates only its block of particles, then all ranks exchange
// (position, velocity, mass) via Allgatherv (7 floats / particle).

void BarnesHutSimulation::integrateMPI_(std::vector<particle_t>& particles,
                                        const SimParams& params)
{
    const float dt = params.dt;

    for (auto& p : particles) {
        p.velocity  += p.acceleration * dt;
        p.position  += p.velocity * dt;
        // acceleration will be overwritten on next computeAccelerationsMPI_
    }
}

#endif // USE_MPI
