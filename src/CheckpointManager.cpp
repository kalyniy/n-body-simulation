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
    std::ofstream file(this->filePath, std::ios::out | std::ios::trunc | std::ios::binary);
    
    if (!file) {
        perror("Error opening file for writing header\n");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        exit(EXIT_FAILURE);
    }

    file.write(reinterpret_cast<const char*>(&header), sizeof(SimulationOutputHeader));
    file.close();
}

void CheckpointManager::increment_passed_steps()
{
    std::fstream file(this->filePath, std::ios::in | std::ios::out | std::ios::binary);
    
    if (!file) {
        perror("Error opening file for incrementing steps\n");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        exit(EXIT_FAILURE);
    }

    // Read the current header
    SimulationOutputHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(SimulationOutputHeader));
    
    // Increment passed_steps
    header.passed_steps++;
    
    // Seek back to the beginning and write updated header
    file.seekp(0, std::ios::beg);
    file.write(reinterpret_cast<const char*>(&header), sizeof(SimulationOutputHeader));
    
    file.close();
}

void CheckpointManager::write_step(particle_t* particles, int count)
{
    std::fstream file(this->filePath, std::ios::in | std::ios::out | std::ios::binary | std::ios::ate);
    
    if (!file) {
        perror("Error opening file for writing step\n");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        exit(EXIT_FAILURE);
    }

    // Flatten particle positions into array
    // Each particle contributes 3 floats (x, y, z)
    size_t flat_size = count * 3;
    float* flat_positions = new float[flat_size];
    
    for (int i = 0; i < count; i++) {
        flat_positions[i * 3 + 0] = particles[i].position.x;
        flat_positions[i * 3 + 1] = particles[i].position.y;
        flat_positions[i * 3 + 2] = particles[i].position.z;
    }
    
    // Write the flattened positions to the end of the file
    file.write(reinterpret_cast<const char*>(flat_positions), flat_size * sizeof(float));
    
    delete[] flat_positions;
    file.close();
    
    // Update the passed_steps counter
    this->increment_passed_steps();
}

SimulationOutputHeader CheckpointManager::read_header()
{
    std::ifstream file(this->filePath, std::ios::in | std::ios::binary);
    
    if (!file) {
        perror("Error opening file for reading header\n");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        exit(EXIT_FAILURE);
    }

    SimulationOutputHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(SimulationOutputHeader));
    file.close();
    return header;
}

size_t CheckpointManager::read_step(float* positions_out, size_t step_index)
{
    std::ifstream file(this->filePath, std::ios::in | std::ios::binary);
    
    if (!file) {
        perror("Error opening file for reading step\n");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        exit(EXIT_FAILURE);
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
    // Offset = header size + (step_index * particles_per_step * 3 floats * sizeof(float))
    size_t step_data_size = header.n_particles * 3 * sizeof(float);
    size_t offset = sizeof(SimulationOutputHeader) + (step_index * step_data_size);
    
    file.seekg(offset, std::ios::beg);
    file.read(reinterpret_cast<char*>(positions_out), step_data_size);
    
    file.close();
    return header.n_particles;
}