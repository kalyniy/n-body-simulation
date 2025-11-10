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
    int processes_count; // how many nodes?
    int rank; // which node am i?

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    SimParams params;
    params.G = 1.0f;
    params.dt = 0.5f;        // Smaller timestep for stability!
    params.min_r2 = 2.0f;     // Slightly larger softening

    NBodySimulation sim = NBodySimulation(std::make_unique<NaiveSimulation>(), params);

    // Basic CLI: headless_main [--hacc DIR] [--solar] [--random N] [--steps S] [--out file]
    std::string out = "output.txt";
    std::string hacc_dir;
    std::size_t steps = 1000, randomN = 0;
    bool use_solar = true;

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
        }
        else if (a == "--G" && i + 1 < argc)
        {
            params.G = std::stof(argv[++i]);
        }
    }

    if (!hacc_dir.empty())
    {
        std::cout << "Loading HACC snapshot: " << hacc_dir << "\n";
        sim.loadHACC(hacc_dir);
    }
    else if (randomN > 0)
    {
        const int nitems = 4;
        int blocklengths[4] = {3, 3, 3, 1};  // 3 floats for each vector3_t, 1 for mass
        
        MPI_Datatype types[4] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
        MPI_Datatype mpi_particle_type;
        MPI_Aint offsets[4];

        offsets[0] = offsetof(particle_t, position);
        offsets[1] = offsetof(particle_t, velocity);
        offsets[2] = offsetof(particle_t, acceleration);
        offsets[3] = offsetof(particle_t, mass);
        
        MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_particle_type);
        MPI_Type_commit(&mpi_particle_type);

        if (rank == 0)
        {
            sim.generateRandom(randomN, WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);
            auto particles = sim.particles();
            auto data = particles.data();
            printf("Rank %d: generated particles\n", rank);
        }

        particle_t* data = new particle_t[randomN];

        MPI_Bcast((void*)data, randomN, mpi_particle_type, 0, MPI_COMM_WORLD);
        printf("Rank %d received data after bcast.\n", rank);

        //if (rank != 0)
        {
            particle_t particle = data[0];

            printf("particle[0] (mass: %f)= pos.x: %f, pos.y: %f, pos.z: %f\n", particle.mass, particle.position.x, particle.position.y, particle.position.z);
            fflush(stdout);
        }
    }
    else if (use_solar)
    {
        std::cout << "Setting up solar system\n";
        sim.setupSolarSystem(WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);
    }

    PerformanceLogger logger(out, RunMode::Run);
    if (!logger.ok())
        std::cerr << "Warning: could not open output log: " << out << "\n";

    size_t log_step_size = 100;

    const size_t n_particles = randomN; //sizes[i];
    std::cout << "Generating " << n_particles << " random particles\n";

    //sim.generateRandom(n_particles, WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);

    CheckpointManager* checkpointManager;

    if (rank == 0)
    {
        auto instance = checkpointManager->getInstance();

        instance->setFilePath("data.bin");
        //instance->write_header();
    }


    auto simulation_start = std::chrono::high_resolution_clock::now();

    for (size_t step = 0; step < steps; step++)
    {
        // auto t0 = std::chrono::high_resolution_clock::now();
        sim.step();
        // auto t1 = std::chrono::high_resolution_clock::now();
        // double dt = std::chrono::duration<double>(t1 - t0).count();

        if ((step % log_step_size) == 0)
            std::cout << "Step " << step << "\n";
    }

    auto simulation_end = std::chrono::high_resolution_clock::now();

    double simulation_time = std::chrono::duration<double>(simulation_end - simulation_start).count();

    logger.log(n_particles, steps, simulation_time);
    
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}