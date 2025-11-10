#include "BarnesHutSimulation.h"
#include <cmath>

BarnesHutSimulation::BarnesHutSimulation(float theta, int leafBucketSize, int maxDepth)
    : theta_(theta)
{
    bp_.bucket_size = leafBucketSize;
    bp_.max_depth   = maxDepth;
    bp_.bounds_pad  = 1e-2f;
}

void BarnesHutSimulation::computeStep(std::vector<particle_t>& particles, const SimParams& params)
{
#ifdef USE_MPI
    computeAccelerationsMPI_(particles, params);
    integrateMPI_(particles, params);
#else
    computeAccelerations_(particles, params);
    integrate_(particles, params);
#endif
}

void BarnesHutSimulation::computeAccelerations_(std::vector<particle_t>& particles, const SimParams& params)
{
    for (auto& p : particles) p.acceleration = {0, 0, 0};
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

void BarnesHutSimulation::integrate_(std::vector<particle_t>& particles, const SimParams& params)
{
    const float dt = params.dt;
    for (auto& p : particles) {
        p.velocity += p.acceleration * dt;
        p.position += p.velocity * dt;
    }
}

#ifdef USE_MPI

void BarnesHutSimulation::computeAccelerationsMPI_(std::vector<particle_t>& particles, 
                                                   const SimParams& params)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    const size_t n = particles.size();
    const float G = params.G;
    const float eps2 = params.min_r2;
    
    // Reset accelerations
    for (auto& p : particles) p.acceleration = {0, 0, 0};
    if (particles.empty()) return;
    
    // Determine local particle range
    size_t particles_per_proc = n / size;
    size_t remainder = n % size;
    size_t start = rank * particles_per_proc + std::min((size_t)rank, remainder);
    size_t end = start + particles_per_proc + (rank < remainder ? 1 : 0);
    
    // CRITICAL: Only rank 0 builds the tree
    Octree tree;
    if (rank == 0) {
        tree.build(particles, bp_);
    }
    
    // Broadcast tree is complex, so simpler approach:
    // Rank 0 computes ALL accelerations, then broadcasts them
    if (rank == 0) {
        for (size_t i = 0; i < n; ++i) {
            vector3_t ai = tree.accelerationOn(i, G, eps2, theta_);
            particles[i].acceleration = ai;
        }
    }
    
    // Broadcast accelerations to all processes
    std::vector<float> acc_data(n * 3);
    if (rank == 0) {
        for (size_t i = 0; i < n; ++i) {
            acc_data[i * 3 + 0] = particles[i].acceleration.x;
            acc_data[i * 3 + 1] = particles[i].acceleration.y;
            acc_data[i * 3 + 2] = particles[i].acceleration.z;
        }
    }
    
    MPI_Bcast(acc_data.data(), n * 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    // Unpack accelerations
    for (size_t i = 0; i < n; ++i) {
        particles[i].acceleration.x = acc_data[i * 3 + 0];
        particles[i].acceleration.y = acc_data[i * 3 + 1];
        particles[i].acceleration.z = acc_data[i * 3 + 2];
    }
}

void BarnesHutSimulation::integrateMPI_(std::vector<particle_t>& particles, 
                                       const SimParams& params)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    const float dt = params.dt;
    const size_t n = particles.size();
    
    // Determine local particle range
    size_t particles_per_proc = n / size;
    size_t remainder = n % size;
    size_t start = rank * particles_per_proc + std::min((size_t)rank, remainder);
    size_t end = start + particles_per_proc + (rank < remainder ? 1 : 0);
    
    // Update local particles
    for (size_t i = start; i < end; ++i)
    {
        auto& p = particles[i];
        p.velocity += p.acceleration * dt;
        p.position += p.velocity * dt;
    }
    
    // Gather updated particles
    std::vector<int> recvcounts(size);
    std::vector<int> displs(size);
    
    for (int r = 0; r < size; ++r) {
        size_t r_start = r * particles_per_proc + std::min((size_t)r, remainder);
        size_t r_end = r_start + particles_per_proc + (r < remainder ? 1 : 0);
        recvcounts[r] = (r_end - r_start) * 10;
        displs[r] = r_start * 10;
    }
    
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
    
    std::vector<float> recvbuf(n * 10);
    MPI_Allgatherv(sendbuf.data(), local_count * 10, MPI_FLOAT,
                   recvbuf.data(), recvcounts.data(), displs.data(),
                   MPI_FLOAT, MPI_COMM_WORLD);
    
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

#endif // USE_MPI