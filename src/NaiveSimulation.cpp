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
    size_t start = rank * particles_per_proc;
    size_t end = (rank == size - 1) ? n : (rank + 1) * particles_per_proc;

    // Compute interactions for local particles
    for (size_t i = start; i < end; ++i)
    {
        particle_t &a = particles[i];

        // Interactions with particles j > i
        for (size_t j = i + 1; j < n; ++j)
        {
            particle_t &b = particles[j];
            vector3_t r = b.position - a.position;
            float r2 = r.x * r.x + r.y * r.y + r.z * r.z + eps2;

            float inv_r = 1.0f / std::sqrt(r2);
            float inv_r3 = inv_r * inv_r * inv_r;

            vector3_t acc = r * (G * inv_r3);
            vector3_t a_acc_d = acc * b.mass;
            vector3_t b_acc_d = acc * a.mass;
            a.acceleration += a_acc_d;
            b.acceleration -= b_acc_d;
        }
    }

    // Flatten accelerations into a buffer for MPI reduction
    std::vector<float> local_acc(n * 3);
    std::vector<float> global_acc(n * 3);
    
    for (size_t i = 0; i < n; ++i) {
        local_acc[i * 3 + 0] = particles[i].acceleration.x;
        local_acc[i * 3 + 1] = particles[i].acceleration.y;
        local_acc[i * 3 + 2] = particles[i].acceleration.z;
    }

    // Reduce accelerations across all processes
    MPI_Allreduce(local_acc.data(), global_acc.data(), 
                  n * 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    // Update particles with global accelerations
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

    // Determine local particle range
    size_t n = particles.size();
    size_t particles_per_proc = n / size;
    size_t start = rank * particles_per_proc;
    size_t end = (rank == size - 1) ? n : (rank + 1) * particles_per_proc;

    // Update local particles
    for (size_t i = start; i < end; ++i)
    {
        auto &p = particles[i];
        p.velocity += p.acceleration * dt;
        p.position += p.velocity * dt;
    }

    // Flatten particle data for broadcasting
    std::vector<float> particle_data(n * 10); // 3 pos + 3 vel + 3 acc + 1 mass = 10 floats per particle
    
    for (size_t i = start; i < end; ++i) {
        size_t idx = i * 10;
        particle_data[idx + 0] = particles[i].position.x;
        particle_data[idx + 1] = particles[i].position.y;
        particle_data[idx + 2] = particles[i].position.z;
        particle_data[idx + 3] = particles[i].velocity.x;
        particle_data[idx + 4] = particles[i].velocity.y;
        particle_data[idx + 5] = particles[i].velocity.z;
        particle_data[idx + 6] = particles[i].acceleration.x;
        particle_data[idx + 7] = particles[i].acceleration.y;
        particle_data[idx + 8] = particles[i].acceleration.z;
        particle_data[idx + 9] = particles[i].mass;
    }

    // Gather all particle data
    std::vector<float> all_particle_data(n * 10);
    MPI_Allreduce(particle_data.data(), all_particle_data.data(),
                  n * 10, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    // Unpack particle data
    for (size_t i = 0; i < n; ++i) {
        size_t idx = i * 10;
        particles[i].position.x = all_particle_data[idx + 0];
        particles[i].position.y = all_particle_data[idx + 1];
        particles[i].position.z = all_particle_data[idx + 2];
        particles[i].velocity.x = all_particle_data[idx + 3];
        particles[i].velocity.y = all_particle_data[idx + 4];
        particles[i].velocity.z = all_particle_data[idx + 5];
        particles[i].acceleration.x = all_particle_data[idx + 6];
        particles[i].acceleration.y = all_particle_data[idx + 7];
        particles[i].acceleration.z = all_particle_data[idx + 8];
        particles[i].mass = all_particle_data[idx + 9];
    }
}