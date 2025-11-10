#include "CheckpointManager.h"
#include <fstream>
#ifdef USE_MPI
#include <mpi.h>    
#endif

void CheckpointManager::write_header(SimulationOutputHeader header)
{
    std::fstream file(this->filePath, std::ios::in | std::ios::out | std::ios::trunc | std::ios::binary);
    
    if (!file) {
        perror("Error opening file\n");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        exit(EXIT_FAILURE);
    }

    file.write((char*) &header, sizeof(SimulationOutputHeader));
    file.close();
}

void CheckpointManager::write_step(particle_t* particles, int count)
{
    this->increment_passed_steps();
}

SimulationOutputHeader CheckpointManager::read_header()
{
    std::fstream file(this->filePath, std::ios::in | std::ios::out | std::ios::binary);
    
    if (!file) {
        perror("Error opening file\n");
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        exit(EXIT_FAILURE);
    }

    SimulationOutputHeader header;
    file.read((char*) &header, sizeof(SimulationOutputHeader));
    file.close();
    return header;
}