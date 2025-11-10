#include "NaiveSimulation.h"

#include <stdlib.h>
#include <cmath>
#include <vector>

#include <mpi.h>

void NaiveSimulation::computeStep(std::vector<particle_t>& particles, const SimParams& params) {
    computeAccelerations_(particles, params);
    integrate_(particles, params);
}

void NaiveSimulation::computeAccelerations_(std::vector<particle_t>& particles, const SimParams& params)
{
    const float eps2 = params.min_r2;
    const float G = params.G;
    const size_t n = particles.size();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Clear accelerations
    for (auto &p : particles)
        p.acceleration = {0, 0, 0};

    // Determine local particle range for this process
    size_t particles_per_proc = n / size;
    size_t remainder = n % size;
    size_t start = rank * particles_per_proc + std::min((size_t)rank, remainder);
    size_t end = start + particles_per_proc + (rank < remainder ? 1 : 0);

    // Compute interactions for local particles
    for (size_t i = start; i < end; ++i)
    {
        particle_t &a = particles[i];

        // Compute all interactions for particle i
        for (size_t j = 0; j < n; ++j)
        {
            if (i == j) continue;
            
            particle_t &b = particles[j];
            vector3_t r = b.position - a.position;
            float r2 = r.x * r.x + r.y * r.y + r.z * r.z + eps2;

            float inv_r = 1.0f / std::sqrt(r2);
            float inv_r3 = inv_r * inv_r * inv_r;

            vector3_t acc = r * (G * inv_r3 * b.mass);
            a.acceleration += acc;
        }
    }

    // Gather accelerations from all processes
    std::vector<float> local_acc(n * 3, 0.0f);
    std::vector<float> global_acc(n * 3);
    
    // Pack only local accelerations
    for (size_t i = start; i < end; ++i) {
        local_acc[i * 3 + 0] = particles[i].acceleration.x;
        local_acc[i * 3 + 1] = particles[i].acceleration.y;
        local_acc[i * 3 + 2] = particles[i].acceleration.z;
    }

    // Reduce accelerations (sum contributions from all processes)
    MPI_Allreduce(local_acc.data(), global_acc.data(), 
                  n * 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    // Update all particles with global accelerations
    for (size_t i = 0; i < n; ++i) {
        particles[i].acceleration.x = global_acc[i * 3 + 0];
        particles[i].acceleration.y = global_acc[i * 3 + 1];
        particles[i].acceleration.z = global_acc[i * 3 + 2];
    }
}

void NaiveSimulation::integrate_(std::vector<particle_t>& particles, const SimParams& params)
{
    const float dt = params.dt;
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t n = particles.size();
    size_t particles_per_proc = n / size;
    size_t remainder = n % size;
    size_t start = rank * particles_per_proc + std::min((size_t)rank, remainder);
    size_t end = start + particles_per_proc + (rank < remainder ? 1 : 0);

    // Update local particles
    for (size_t i = start; i < end; ++i)
    {
        auto &p = particles[i];
        p.velocity += p.acceleration * dt;
        p.position += p.velocity * dt;
    }

    // Gather updated particles from all processes
    // Use Allgatherv to handle uneven distribution
    std::vector<int> recvcounts(size);
    std::vector<int> displs(size);
    
    for (int r = 0; r < size; ++r) {
        size_t r_start = r * particles_per_proc + std::min((size_t)r, remainder);
        size_t r_end = r_start + particles_per_proc + (r < remainder ? 1 : 0);
        recvcounts[r] = (r_end - r_start) * 10; // 10 floats per particle
        displs[r] = r_start * 10;
    }

    // Pack local particle data
    size_t local_count = end - start;
    std::vector<float> sendbuf(local_count * 10);
    
    for (size_t i = 0; i < local_count; ++i) {
        size_t pi = start + i;
        size_t idx = i * 10;
        sendbuf[idx + 0] = particles[pi].position.x;
        sendbuf[idx + 1] = particles[pi].position.y;
        sendbuf[idx + 2] = particles[pi].position.z;
        sendbuf[idx + 3] = particles[pi].velocity.x;
        sendbuf[idx + 4] = particles[pi].velocity.y;
        sendbuf[idx + 5] = particles[pi].velocity.z;
        sendbuf[idx + 6] = particles[pi].acceleration.x;
        sendbuf[idx + 7] = particles[pi].acceleration.y;
        sendbuf[idx + 8] = particles[pi].acceleration.z;
        sendbuf[idx + 9] = particles[pi].mass;
    }

    // Gather all particle data
    std::vector<float> recvbuf(n * 10);
    MPI_Allgatherv(sendbuf.data(), local_count * 10, MPI_FLOAT,
                   recvbuf.data(), recvcounts.data(), displs.data(),
                   MPI_FLOAT, MPI_COMM_WORLD);

    // Unpack all particle data
    for (size_t i = 0; i < n; ++i) {
        size_t idx = i * 10;
        particles[i].position.x = recvbuf[idx + 0];
        particles[i].position.y = recvbuf[idx + 1];
        particles[i].position.z = recvbuf[idx + 2];
        particles[i].velocity.x = recvbuf[idx + 3];
        particles[i].velocity.y = recvbuf[idx + 4];
        particles[i].velocity.z = recvbuf[idx + 5];
        particles[i].acceleration.x = recvbuf[idx + 6];
        particles[i].acceleration.y = recvbuf[idx + 7];
        particles[i].acceleration.z = recvbuf[idx + 8];
        particles[i].mass = recvbuf[idx + 9];
    }
}