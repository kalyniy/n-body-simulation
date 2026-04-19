#include "CheckpointManager.h"
#include <fstream>
#ifdef USE_MPI
#include <mpi.h>    
#endif

// Initialize static members
CheckpointManager* CheckpointManager::instancePtr = nullptr;
std::mutex CheckpointManager::mtx;

void CheckpointManager::write_header(SimulationOutputHeader header)
{
    std::ofstream file(filePath, std::ios::out | std::ios::trunc | std::ios::binary);
    
    if (!file) {
        perror("Error writing header");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        std::exit(EXIT_FAILURE);
    }

    file.write(reinterpret_cast<const char*>(&header), sizeof(SimulationOutputHeader));
    file.close();
}

void CheckpointManager::write_masses(const ParticleSoA& particles)
{
    std::fstream file(filePath, std::ios::in | std::ios::out | std::ios::binary | std::ios::ate);
    
    if (!file) { 
        perror("Error writing masses");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        std::exit(EXIT_FAILURE);
    }

    file.write(reinterpret_cast<const char*>(particles.mass.data()), 
               particles.size() * sizeof(float));
    
    file.close();
}

void CheckpointManager::increment_passed_steps()
{
    std::fstream file(filePath, std::ios::in | std::ios::out | std::ios::binary);
    if (!file) {
        perror("Error incrementing steps");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        std::exit(EXIT_FAILURE);
    }

    SimulationOutputHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(SimulationOutputHeader));
    
    // Increment passed_steps
    header.passed_steps++;
    
    // Seek back to the beginning and write updated header
    file.seekp(0, std::ios::beg);
    file.write(reinterpret_cast<const char*>(&header), sizeof(SimulationOutputHeader));
    
    file.close();
}

void CheckpointManager::write_step(const ParticleSoA& particles)
{
    std::fstream file(filePath, std::ios::in | std::ios::out | std::ios::binary | std::ios::ate);
    if (!file) {
        perror("Error writing step");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        std::exit(EXIT_FAILURE);
    }

    size_t n = particles.size();
    float pos[3] = { 0.0f, 0.0f, 0.0f };

    // Write interleaved x,y,z positions
    for (size_t i = 0; i < n; ++i) {
        pos[0] = particles.pos_x[i];
        pos[1] = particles.pos_y[i];
        pos[2] = particles.pos_z[i];
        file.write(reinterpret_cast<const char*>(pos), sizeof(pos));
    }

    file.close();

    increment_passed_steps();
}

SimulationOutputHeader CheckpointManager::read_header()
{
    std::ifstream file(filePath, std::ios::in | std::ios::binary);

    if (!file) { 
        perror("Error reading header");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        std::exit(EXIT_FAILURE);
    }

    SimulationOutputHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(SimulationOutputHeader));
    file.close();
    return header;
}

void CheckpointManager::read_masses(float* masses_out, size_t n_particles)
{
    std::ifstream file(filePath, std::ios::in | std::ios::binary);
    if (!file) {
        perror("Error reading masses");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        std::exit(EXIT_FAILURE);
    }

    file.seekg(sizeof(SimulationOutputHeader), std::ios::beg);
    file.read(reinterpret_cast<char*>(masses_out), n_particles * sizeof(float));
    file.close();
}

size_t CheckpointManager::read_step(float* positions_out, size_t step_index)
{
    std::ifstream file(filePath, std::ios::in | std::ios::binary);
    if (!file) {
        perror("Error reading step");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        std::exit(EXIT_FAILURE);
    }

    // Read header to get particle count
    SimulationOutputHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(SimulationOutputHeader));

    // Check if step_index is valid
    if (step_index >= header.passed_steps) {
        file.close();
        return 0; // No data available for this step yet
    }
    
    // Calculate offset to the desired step
    // Offset = header + masses + (step_index * positions_per_step)
    size_t mass_data_size = header.n_particles * sizeof(float);
    size_t step_data_size = header.n_particles * 3 * sizeof(float);
    size_t offset = sizeof(SimulationOutputHeader) + mass_data_size + (step_index * step_data_size);
    
    file.seekg(offset, std::ios::beg);
    file.read(reinterpret_cast<char*>(positions_out), step_data_size);
    
    file.close();
    return header.n_particles;
}