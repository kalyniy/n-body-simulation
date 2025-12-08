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

    const size_t n    = particles.size();
    const float  G    = params.G;
    const float  eps2 = params.min_r2;

    // Reset accelerations
    for (auto& p : particles) {
        p.acceleration = {0.0f, 0.0f, 0.0f};
    }
    if (n == 0) return;

    int nodeCount      = 0;
    int leafIndexCount = 0;

    // 1) Rank 0 builds the Octree and exports to MPI-friendly arrays
    if (rank == 0) {
        Octree tree;
        tree.build(particles, bp_);
        tree.exportToMpiTree(mpi_nodes_, mpi_leafIndices_);

        nodeCount      = static_cast<int>(mpi_nodes_.size());
        leafIndexCount = static_cast<int>(mpi_leafIndices_.size());
    }

    // 2) Broadcast tree sizes
    MPI_Bcast(&nodeCount,      1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&leafIndexCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    mpi_nodes_.resize(nodeCount);
    mpi_leafIndices_.resize(leafIndexCount);

    // 3) Broadcast tree data
    if (nodeCount > 0) {
        MPI_Bcast(reinterpret_cast<unsigned char*>(mpi_nodes_.data()),
                  nodeCount * static_cast<int>(sizeof(Octree::MpiTreeNode)),
                  MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    if (leafIndexCount > 0) {
        MPI_Bcast(mpi_leafIndices_.data(),
                  leafIndexCount,
                  MPI_INT, 0, MPI_COMM_WORLD);
    }

    // 4) Determine local particle range for this rank
    const size_t base      = n / size;
    const size_t remainder = n % size;

    const size_t start =
        static_cast<size_t>(rank) * base +
        std::min(static_cast<size_t>(rank), remainder);
    const size_t end =
        start + base + (rank < static_cast<int>(remainder) ? 1 : 0);

    // 5) Compute accelerations for local subset
    for (size_t i = start; i < end; ++i) {
        vector3_t acc = {0.0f, 0.0f, 0.0f};
        if (nodeCount > 0) {
            traverseMpiTreeAccumulate_(0, static_cast<int>(i),
                                       G, eps2, theta_,
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
                int j = mpi_leafIndices_[offset + k];
                if (j == i) continue;

                const particle_t& pj = particles[j];
                vector3_t rij = {
                    pj.position.x - pi.x,
                    pj.position.y - pi.y,
                    pj.position.z - pi.z
                };
                float d2     = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z + eps2;
                float inv_r  = 1.0f / std::sqrt(d2);
                float inv_r3 = inv_r * inv_r * inv_r;
                float s      = G * pj.mass * inv_r3;
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
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const float  dt = params.dt;
    const size_t n  = particles.size();
    if (n == 0) return;

    const size_t base      = n / size;
    const size_t remainder = n % size;

    const size_t start =
        static_cast<size_t>(rank) * base +
        std::min(static_cast<size_t>(rank), remainder);
    const size_t end =
        start + base + (rank < static_cast<int>(remainder) ? 1 : 0);
    const size_t local_count = (end > start) ? (end - start) : 0;

    // 1) Update our local chunk
    for (size_t i = start; i < end; ++i) {
        auto& p = particles[i];
        p.velocity  += p.acceleration * dt;
        p.position  += p.velocity * dt;
    }

    // 2) Pack local data: pos (3) + vel (3) + mass (1) = 7 floats
    std::vector<float> sendbuf(local_count * 7);
    for (size_t i = 0; i < local_count; ++i) {
        const auto& p = particles[start + i];
        const size_t idx = i * 7;
        sendbuf[idx + 0] = p.position.x;
        sendbuf[idx + 1] = p.position.y;
        sendbuf[idx + 2] = p.position.z;
        sendbuf[idx + 3] = p.velocity.x;
        sendbuf[idx + 4] = p.velocity.y;
        sendbuf[idx + 5] = p.velocity.z;
        sendbuf[idx + 6] = p.mass;
    }

    // 3) Prepare recvcounts / displs for Allgatherv
    std::vector<int> recvcounts(size);
    std::vector<int> displs(size);

    for (int r = 0; r < size; ++r) {
        const size_t r_base      = n / size;
        const size_t r_remainder = n % size;
        const size_t r_start =
            static_cast<size_t>(r) * r_base +
            std::min(static_cast<size_t>(r), r_remainder);
        const size_t r_end =
            r_start + r_base +
            (r < static_cast<int>(r_remainder) ? 1 : 0);
        const size_t r_local_count = (r_end > r_start) ? (r_end - r_start) : 0;

        recvcounts[r] = static_cast<int>(r_local_count * 7);
        displs[r]     = static_cast<int>(r_start * 7);
    }

    // 4) Gather updated state to all ranks
    std::vector<float> recvbuf(n * 7);
    MPI_Allgatherv(sendbuf.data(),
                   static_cast<int>(sendbuf.size()),
                   MPI_FLOAT,
                   recvbuf.data(),
                   recvcounts.data(),
                   displs.data(),
                   MPI_FLOAT,
                   MPI_COMM_WORLD);

    // 5) Unpack into our local particle array (replicated on all ranks)
    for (size_t i = 0; i < n; ++i) {
        const size_t idx = i * 7;
        auto& p = particles[i];
        p.position.x = recvbuf[idx + 0];
        p.position.y = recvbuf[idx + 1];
        p.position.z = recvbuf[idx + 2];
        p.velocity.x = recvbuf[idx + 3];
        p.velocity.y = recvbuf[idx + 4];
        p.velocity.z = recvbuf[idx + 5];
        p.mass       = recvbuf[idx + 6];
        // acceleration will be overwritten next step
    }
}

#endif // USE_MPI
