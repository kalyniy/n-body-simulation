#include "MpiNaiveSimulation.h"
#include <cmath>
#include <vector>
#include <algorithm>
#ifdef USE_MPI
#include <mpi.h>
#endif
void MpiNaiveSimulation::computeStep(ParticleSoA& particles, const SimParams& params) {
    computeAccelerations_(particles, params);
    integrate_(particles, params);
}

void MpiNaiveSimulation::computeAccelerations_(ParticleSoA& p, const SimParams& params)
{
#ifdef USE_MPI
    const float eps2 = params.min_r2;
    const float G    = params.G;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int n_local = static_cast<int>(p.size());
    p.zeroAccelerations();

    // Gather counts
    std::vector<int> counts(size);
    MPI_Allgather(&n_local, 1, MPI_INT,
                  counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    int max_n = *std::max_element(counts.begin(), counts.end());

    // Ring buffers: 4 floats per particle (px, py, pz, mass)
    std::vector<float> B(4 * max_n);     // current block
    std::vector<float> B_in(4 * max_n);  // receive buffer

    // Pack our own particles
    for (int i = 0; i < n_local; ++i) {
        B[4*i + 0] = p.pos_x[i];
        B[4*i + 1] = p.pos_y[i];
        B[4*i + 2] = p.pos_z[i];
        B[4*i + 3] = p.mass[i];
    }

    int prev = (rank - 1 + size) % size;
    int next = (rank + 1) % size;

    // 4) Ring over all ranks
    for (int step = 0; step < size; ++step) {
        // Who owns the block currently in B?
        // At step 0: owner = rank
        // At step 1: owner = (rank - 1 + size) % size
        // ...
        int owner = (rank - step + size) % size;
        int src_count = counts[owner]; // number of valid particles in B

        for (int i = 0; i < n_local; ++i) {
            for (int j = 0; j < src_count; ++j) {
                if (owner == rank && i == j) continue;

                float rx = B[4*j + 0] - p.pos_x[i];
                float ry = B[4*j + 1] - p.pos_y[i];
                float rz = B[4*j + 2] - p.pos_z[i];

                float r2 = rx * rx + ry * ry + rz * rz + eps2;
                float inv_r  = 1.0f / std::sqrt(r2);
                float inv_r3 = inv_r * inv_r * inv_r;
                float s = G * B[4*j + 3] * inv_r3;

                p.acc_x[i] += s * rx;
                p.acc_y[i] += s * ry;
                p.acc_z[i] += s * rz;
            }
        }
        if (step < size - 1) {
            int next_owner = (owner - 1 + size) % size;
            int recv_count = counts[next_owner];

            MPI_Sendrecv(
                B.data(),    4 * src_count,  MPI_FLOAT, next, 0,
                B_in.data(), 4 * recv_count, MPI_FLOAT, prev, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            B.swap(B_in);
        }
    }
#else
    (void)p;
    (void)params;
    std::exit(EXIT_FAILURE);
#endif
}

void MpiNaiveSimulation::integrate_(ParticleSoA& p, const SimParams& params)
{
    const float dt = params.dt;
    for (size_t i = 0; i < p.size(); ++i) {
        p.vel_x[i] += p.acc_x[i] * dt;
        p.vel_y[i] += p.acc_y[i] * dt;
        p.vel_z[i] += p.acc_z[i] * dt;
        p.pos_x[i] += p.vel_x[i] * dt;
        p.pos_y[i] += p.vel_y[i] * dt;
        p.pos_z[i] += p.vel_z[i] * dt;
    }
}