#include "MpiNaiveSimulation.h"

#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>
#ifdef USE_MPI
#include <mpi.h>
#endif
void MpiNaiveSimulation::computeStep(std::vector<particle_t>& particles, const SimParams& params) {
    computeAccelerations_(particles, params);
    integrate_(particles, params);
}

void MpiNaiveSimulation::computeAccelerations_(std::vector<particle_t>& particles, const SimParams& params)
{
#ifdef USE_MPI
    const float eps2 = params.min_r2;
    const float G    = params.G;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int n_local = static_cast<int>(particles.size());

    // 1) Clear local accelerations
    for (auto &p : particles)
        p.acceleration = {0.0f, 0.0f, 0.0f};

    // 2) Every rank tells others how many particles it owns (local subset)
    std::vector<int> counts(size);
    MPI_Allgather(&n_local, 1, MPI_INT,
                  counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    // max local size (for buffers)
    int max_n_local = *std::max_element(counts.begin(), counts.end());

    // 3) Allocate ring buffers (sources j)
    std::vector<particle_t> B(max_n_local);     // current block
    std::vector<particle_t> B_in(max_n_local);  // receive buffer

    // Copy our own local particles into B initially
    std::copy(particles.begin(), particles.end(), B.begin());

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

        // 4a) Compute contributions from block B onto our local particles
        if (owner == rank) {
            // We are using our own block; need to avoid self-interaction
            for (int i = 0; i < n_local; ++i) {
                particle_t &a = particles[i];
                for (int j = 0; j < n_local; ++j) {
                    if (i == j) continue;  // skip self
                    const particle_t &b = B[j];

                    vector3_t r = { b.position.x - a.position.x,
                                    b.position.y - a.position.y,
                                    b.position.z - a.position.z };

                    float r2 = r.x * r.x + r.y * r.y + r.z * r.z + eps2;
                    float inv_r  = 1.0f / std::sqrt(r2);
                    float inv_r3 = inv_r * inv_r * inv_r;

                    float s = G * b.mass * inv_r3;
                    a.acceleration.x += s * r.x;
                    a.acceleration.y += s * r.y;
                    a.acceleration.z += s * r.z;
                }
            }
        } else {
            // Block belongs to some other rank; no self-interactions to worry about
            for (int i = 0; i < n_local; ++i) {
                particle_t &a = particles[i];
                for (int j = 0; j < src_count; ++j) {
                    const particle_t &b = B[j];

                    vector3_t r = { b.position.x - a.position.x,
                                    b.position.y - a.position.y,
                                    b.position.z - a.position.z };

                    float r2 = r.x * r.x + r.y * r.y + r.z * r.z + eps2;
                    float inv_r  = 1.0f / std::sqrt(r2);
                    float inv_r3 = inv_r * inv_r * inv_r;

                    float s = G * b.mass * inv_r3;
                    a.acceleration.x += s * r.x;
                    a.acceleration.y += s * r.y;
                    a.acceleration.z += s * r.z;
                }
            }
        }

        // 4b) Rotate the block B, except after the last step
        if (step < size - 1) {
            // After this rotation, B will hold the block from owner' = (owner - 1 + size) % size
            int next_owner = (owner - 1 + size) % size;
            int recv_count = counts[next_owner];

            // Send src_count particles, receive recv_count particles
            MPI_Sendrecv(
                B.data(), src_count * static_cast<int>(sizeof(particle_t)), MPI_BYTE, next, 0,
                B_in.data(), recv_count * static_cast<int>(sizeof(particle_t)), MPI_BYTE, prev, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );

            // Swap buffers so B holds the newly received block
            B.swap(B_in);
        }
    }
#else
    exit(EXIT_FAILURE);
#endif
}

void MpiNaiveSimulation::integrate_(std::vector<particle_t>& particles, const SimParams& params)
{
    const float dt = params.dt;

    // No MPI calls needed; each rank updates only its local chunk.
    for (auto &p : particles) {
        p.velocity.x += p.acceleration.x * dt;
        p.velocity.y += p.acceleration.y * dt;
        p.velocity.z += p.acceleration.z * dt;

        p.position.x += p.velocity.x * dt;
        p.position.y += p.velocity.y * dt;
        p.position.z += p.velocity.z * dt;
    }
}