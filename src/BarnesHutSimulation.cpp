// BarnesHutSimulation.cpp
#include "BarnesHutSimulation.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

BarnesHutSimulation::BarnesHutSimulation(float theta,
                                         int   leafBucketSize,
                                         int   maxDepth)
    : theta_(theta)
{
    bp_.bucket_size = leafBucketSize;
    bp_.max_depth   = maxDepth;
    bp_.bounds_pad  = 1e-2f;
}

void BarnesHutSimulation::computeStep(ParticleSoA& particles, const SimParams& params)
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

void BarnesHutSimulation::computeAccelerations_(ParticleSoA& particles, const SimParams& params)
{
    particles.zeroAccelerations();
    if (particles.empty()) return;

    int n = static_cast<int>(p.size());
    Octree tree;
    tree.build(particles.pos_x.data(), particles.pos_y.data(), particles.pos_z.data(),
               particles.mass.data(), n, bp_);

    for (int i = 0; i < n; ++i)
    {
        tree.accelerationOn(i, params.G, params.min_r2, theta_,
                            particles.acc_x[i], particles.acc_y[i], particles.acc_z[i]);
    }
}

void BarnesHutSimulation::integrate_(ParticleSoA& particles, const SimParams& params)
{
    const float dt = params.dt;
    for (size_t i = 0; i < particles.size(); ++i) {
        particles.vel_x[i] += particles.acc_x[i] * dt;
        particles.vel_y[i] += particles.acc_y[i] * dt;
        particles.vel_z[i] += particles.acc_z[i] * dt;
        particles.pos_x[i] += particles.vel_x[i] * dt;
        particles.pos_y[i] += particles.vel_y[i] * dt;
        particles.pos_z[i] += particles.vel_z[i] * dt;
    }
}

#ifdef USE_MPI

// ======================= MPI: Accelerations =======================
//
// Rank 0 builds the tree, flattens it, and broadcasts to all ranks.
// All ranks traverse the same tree but only for their subset of particles.
void BarnesHutSimulation::computeAccelerationsMPI_(ParticleSoA& p,
                                                   const SimParams& params)
{
    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    const int   localN = static_cast<int>(p.size());
    p.zeroAccelerations();
    if (localN == 0) return;

    // Gather counts
    std::vector<int> counts(nproc);
    MPI_Allgather(&localN, 1, MPI_INT,
                  counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    int globalN = 0;
    std::vector<int> displs(nproc, 0);
    if (rank == 0) {
        for (int r = 0; r < nproc; ++r) {
            displs[r] = globalN;
            globalN  += counts[r];
        }
    }

    // Pack local pos+mass (4 floats per particle)
    std::vector<float> sendbuf(4 * localN);
    for (int i = 0; i < localN; ++i) {
        sendbuf[4*i + 0] = p.pos_x[i];
        sendbuf[4*i + 1] = p.pos_y[i];
        sendbuf[4*i + 2] = p.pos_z[i];
        sendbuf[4*i + 3] = p.mass[i];
    }

    // 3) Gatherv packed data on rank 0 (like reference code's MPI_Gather on rr)
    std::vector<int> recvcountsF;
    std::vector<int> displsF;
    if (rank == 0) {
        recvcountsF.resize(nproc);
        displsF.resize(nproc);

        int offsetF = 0;
        for (int r = 0; r < nproc; ++r) {
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

    // Rank 0 builds tree
    int nodeCount = 0;
    int leafCount = 0;
    if (rank == 0) {
        // Unpack into temporary SoA for tree build
        std::vector<float> gpx(globalN), gpy(globalN), gpz(globalN), gm(globalN);
        for (int i = 0; i < globalN; ++i) {
            gpx[i] = globalPosMass_[4*i + 0];
            gpy[i] = globalPosMass_[4*i + 1];
            gpz[i] = globalPosMass_[4*i + 2];
            gm[i]  = globalPosMass_[4*i + 3];
        }
        Octree tree;
        tree.build(gpx.data(), gpy.data(), gpz.data(), gm.data(), globalN, bp_);
        tree.exportToFlatTree(mpi_nodes_, mpi_leaves_);
        nodeCount = static_cast<int>(mpi_nodes_.size());
        leafCount = static_cast<int>(mpi_leaves_.size());
    }

    // Broadcast tree
    MPI_Bcast(&nodeCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&leafCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    mpi_nodes_.resize(nodeCount);
    mpi_leaves_.resize(leafCount);
    if (nodeCount > 0)
        MPI_Bcast(mpi_nodes_.data(), nodeCount * sizeof(Octree::FlatTreeNode), MPI_BYTE, 0, MPI_COMM_WORLD);
    if (leafCount > 0)
        MPI_Bcast(mpi_leaves_.data(), leafCount * sizeof(Octree::FlatLeafParticle), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Each rank traverses for its local particles
    for (int i = 0; i < localN; ++i) {
        float ax = 0, ay = 0, az = 0;
        if (nodeCount > 0)
            traverseMpiTree_(0, i, p.pos_x.data(), p.pos_y.data(), p.pos_z.data(),
                             params.G, params.min_r2, theta_, ax, ay, az);
        p.acc_x[i] = ax; p.acc_y[i] = ay; p.acc_z[i] = az;
    }
}

void BarnesHutSimulation::traverseMpiTree_(
    int nodeId, int i,
    const float* px, const float* py, const float* pz,
    float G, float eps2, float theta,
    float& ax, float& ay, float& az) const
{
    if (nodeId < 0 ||
        nodeId >= static_cast<int>(mpi_nodes_.size()))
    {
        return;
    }

    const auto& node = mpi_nodes_[nodeId];
    if (node.mass <= 0.0f) return;

    float rx = node.com[0] - px[i];
    float ry = node.com[1] - py[i];
    float rz = node.com[2] - pz[i];
    
    float r2 = rx*rx + ry*ry + rz*rz;

    if (node.isLeaf) {
        if (node.leafOffset >= 0) {
            for (int k = 0; k < node.leafCount; ++k) {
                const auto& lp = mpi_leaves_[node.leafOffset + k];
                float dx = lp.pos[0] - px[i], dy = lp.pos[1] - py[i], dz = lp.pos[2] - pz[i];
                float d2 = dx*dx + dy*dy + dz*dz + eps2;
                float inv = 1.f / std::sqrt(d2);
                float inv3 = inv*inv*inv;
                float s = G * lp.mass * inv3;
                ax += dx*s; ay += dy*s; az += dz*s;
            }
        }
        return;
    }

    float s_width = node.half * 2.f;
    float d = std::sqrt(r2) + 1e-12f;
    if ((s_width / d) < theta) {
        float r2s = r2 + eps2;
        float inv = 1.f / std::sqrt(r2s);
        float inv3 = inv* inv* inv;
        float s = G * node.mass * inv3;
        ax += rx*s;
        ay += ry*s;
        az += rz*s;
        return;
    }

    for (int c = 0; c < 8; ++c)
        if (node.child[c] != -1)
            traverseMpiTree_(node.child[c], i, px, py, pz, G, eps2, theta, ax, ay, az);
}

void BarnesHutSimulation::integrateMPI_(ParticleSoA& particles, const SimParams& params)
{
    integrate_(particles, params);
}

#endif // USE_MPI
