// headless_main.cpp
#include <chrono>
#include <iostream>
#include "Simulation.h"
#include "PerformanceLogger.hpp"
#include "NaiveSimulation.h"
#include "BarnesHutSimulation.h"
#include "World.h"
#include "CheckpointManager.h"
#include <mpi.h>

int main(int argc, char **argv)
{
    int processes_count;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    SimParams params;
    params.G = 1.0f;
    params.dt = 0.5f;
    params.min_r2 = 2.0f;

    NBodySimulation sim = NBodySimulation(std::make_unique<NaiveSimulation>(), params);

    std::string out = "output.txt";
    std::string hacc_dir;
    std::size_t steps = 1000, randomN = 0;
    bool use_solar = false;

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (a == "--hacc" && i + 1 < argc)
        {
            hacc_dir = argv[++i];
            use_solar = false;
        }
        else if (a == "--random" && i + 1 < argc)
        {
            randomN = std::stoul(argv[++i]);
            use_solar = false;
        }
        else if (a == "--solar")
        {
            use_solar = true;
        }
        else if (a == "--steps" && i + 1 < argc)
        {
            steps = std::stoul(argv[++i]);
        }
        else if (a == "--out" && i + 1 < argc)
        {
            out = argv[++i];
        }
        else if (a == "--dt" && i + 1 < argc)
        {
            params.dt = std::stof(argv[++i]);
            sim.setDt(params.dt);
        }
        else if (a == "--G" && i + 1 < argc)
        {
            params.G = std::stof(argv[++i]);
            sim.setG(params.G);
        }
    }

    // Create MPI datatype for particle_t
    const int nitems = 4;
    int blocklengths[4] = {3, 3, 3, 1};
    MPI_Datatype types[4] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_particle_type;
    MPI_Aint offsets[4];

    offsets[0] = offsetof(particle_t, position);
    offsets[1] = offsetof(particle_t, velocity);
    offsets[2] = offsetof(particle_t, acceleration);
    offsets[3] = offsetof(particle_t, mass);
    
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_particle_type);
    MPI_Type_commit(&mpi_particle_type);

    size_t n_particles = 0;

    if (!hacc_dir.empty())
    {
        if (rank == 0) {
            std::cout << "Loading HACC snapshot: " << hacc_dir << "\n";
            sim.loadHACC(hacc_dir);
            n_particles = sim.particles().size();
        }
        
        // Broadcast particle count
        MPI_Bcast(&n_particles, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        
        // Prepare buffer
        if (rank != 0) {
            sim.particles().resize(n_particles);
        }
        
        // Broadcast particles
        MPI_Bcast(sim.particles().data(), n_particles, mpi_particle_type, 0, MPI_COMM_WORLD);
    }
    else if (randomN > 0)
    {
        n_particles = randomN;
        
        if (rank == 0)
        {
            //sim.generateRandom(randomN, WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);
            sim.generateGalaxyDisk(n_particles, 
                300.0f,   // disk radius
                50.0f     // disk thickness
            );
            printf("Rank %d: generated %zu particles\n", rank, randomN);
        }
        else
        {
            sim.particles().resize(randomN);
        }

        // Broadcast particles to all processes
        MPI_Bcast(sim.particles().data(), randomN, mpi_particle_type, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            printf("Rank %d: broadcasted particles\n", rank);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        printf("Rank %d: received %zu particles\n", rank, sim.particles().size());
    }
    else if (use_solar)
    {
        if (rank == 0) {
            std::cout << "Setting up solar system\n";
            sim.setupSolarSystem(WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);
            n_particles = sim.particles().size();
        }
        
        // Broadcast particle count
        MPI_Bcast(&n_particles, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        
        // Prepare buffer
        if (rank != 0) {
            sim.particles().resize(n_particles);
        }
        
        // Broadcast particles
        MPI_Bcast(sim.particles().data(), n_particles, mpi_particle_type, 0, MPI_COMM_WORLD);
    }

    PerformanceLogger logger(out, RunMode::Run);
    if (rank == 0 && !logger.ok())
        std::cerr << "Warning: could not open output log: " << out << "\n";

    size_t log_step_size = 100;

    // Only rank 0 handles checkpointing
    if (rank == 0) 
    {
        CheckpointManager* checkpoint = CheckpointManager::getInstance();
        checkpoint->setFilePath("simulation_output.bin");
        
        SimulationOutputHeader header;
        header.n_particles = n_particles;
        header.target_steps = steps;
        header.passed_steps = 0;
        checkpoint->write_header(header);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto simulation_start = std::chrono::high_resolution_clock::now();

    for (size_t step = 0; step < steps; step++)
    {
        sim.step();
        
        if (rank == 0) {
            CheckpointManager* checkpoint = CheckpointManager::getInstance();
            checkpoint->write_step(sim.particles().data(), n_particles);
            
            if (step % log_step_size == 0) {
                std::cout << "Step " << step << " / " << steps << " completed" << std::endl;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto simulation_end = std::chrono::high_resolution_clock::now();

    if (rank == 0) {
        double simulation_time = std::chrono::duration<double>(simulation_end - simulation_start).count();
        logger.log(n_particles, steps, simulation_time);
        std::cout << "Simulation completed in " << simulation_time << " seconds\n";
    }

    MPI_Type_free(&mpi_particle_type);
    MPI_Finalize();
    return 0;
}