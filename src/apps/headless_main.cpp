// headless_main.cpp
#include <chrono>
#include <iostream>
#include <vector>
#include <cstddef> // offsetof

#include "Simulation.h"
#include "PerformanceLogger.hpp"
#include "NaiveSimulation.h"
#include "BarnesHutSimulation.h"
#include "World.h"
#include "CheckpointManager.h"

#ifdef USE_MPI
#include <mpi.h>
#include "MpiNaiveSimulation.h"
#endif

int main(int argc, char **argv)
{
#ifdef USE_MPI
    int processes_count;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);
#else
    int rank = 0;              // for conditional prints
    int processes_count = 1;   // dummy
#endif

    SimParams params;
    params.G      = 1.0f;
    params.dt     = 0.5f;
    params.min_r2 = 2.0f;

    std::unique_ptr<SimulationAlgorithm> algorithm;
    AlgorithmKind algo_kind;
    bool isNaive = false;
    
#ifdef USE_MPI
    if (isNaive)
    {
        std::cout << "MpiNaiveSimulation\n";
        algorithm = std::make_unique<MpiNaiveSimulation>();
        algo_kind  = AlgorithmKind::NaiveMpi;
    }
    else 
    {
        std::cout << "BarnesHutSimulation with MPI\n";
        algorithm = std::make_unique<BarnesHutSimulation>();
        algo_kind  = AlgorithmKind::BarnesHutMpi;
    }
    
#else
    if (isNaive)
    {
        std::cout << "NaiveSimulation\n";
        algorithm = std::make_unique<NaiveSimulation>();
        algo_kind  = AlgorithmKind::NaiveSeq;
    }
    else
    {
        std::cout << "BarnesHutSimulation\n";
        algorithm = std::make_unique<BarnesHutSimulation>();
        algo_kind  = AlgorithmKind::BarnesHutSeq;
    }
#endif

    NBodySimulation sim(std::move(algorithm), params);

    std::string out      = "output.txt";
    std::string hacc_dir;
    std::size_t steps    = 1000;
    std::size_t randomN  = 0;
    bool use_solar       = false;
    bool disable_output = false;

    // -----------------------
    // Parse command line args
    // -----------------------
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
        else if (a == "--disable-output")
        {
            disable_output = true;
        }
    }

#ifdef USE_MPI
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
#endif

    // -------------------------------------------------
    // Global particle count & domain decomposition data
    // -------------------------------------------------
    size_t n_particles_global = 0;

#ifdef USE_MPI
    std::vector<particle_t> global_particles; // only rank 0 uses this
    std::vector<int> counts;
    std::vector<int> displs;
    int n_local = 0;

    // Both NaiveMpi and BarnesHutMpi use index-based domain decomposition
    bool domain_decomposed =
        (algo_kind == AlgorithmKind::NaiveMpi ||
         algo_kind == AlgorithmKind::BarnesHutMpi);
#endif

    // -------------------------
    // Load / generate particles
    // -------------------------
    if (!hacc_dir.empty())
    {
#ifdef USE_MPI
        if (rank == 0)
        {
            std::cout << "Loading HACC snapshot: " << hacc_dir << "\n";
            sim.loadHACC(hacc_dir);
            n_particles_global = sim.particles().size();
            global_particles   = sim.particles(); // copy global array
        }
#else
        std::cout << "Loading HACC snapshot: " << hacc_dir << "\n";
        sim.loadHACC(hacc_dir);
        n_particles_global = sim.particles().size();
#endif
    }
    else if (randomN > 0)
    {
#ifdef USE_MPI
        if (rank == 0)
        {
            n_particles_global = randomN;
            /*
            sim.generateGalaxyDisk(
                n_particles_global,
                300.0f,
                50.0f
            );
            */

            sim.generatePlummerSphere(
                n_particles_global,
                100.0f,     // scale_radius
                (double)n_particles_global * (double)(1.0f)    // total_mass (1.0 per particle)
            );

            std::printf("Rank %d: generated %zu particles\n", rank, n_particles_global);
        }
        #else
        n_particles_global = randomN;
        sim.generateGalaxyDisk(
            n_particles_global,
            300.0f,
            50.0f
        );
#endif
    }
    else if (use_solar)
    {
#ifdef USE_MPI
        if (rank == 0)
        {
            std::cout << "Setting up solar system\n";
            sim.setupSolarSystem(WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);
            n_particles_global = sim.particles().size();
            global_particles   = sim.particles();
        }
#else
        std::cout << "Setting up solar system\n";
        sim.setupSolarSystem(WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);
        n_particles_global = sim.particles().size();
#endif
    }
    else
    {
        // If nothing specified, you might want to default to something
        if (rank == 0)
            std::cerr << "No initial condition specified (--hacc, --random, or --solar)\n";
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#else
        return 1;
#endif
    }

#ifdef USE_MPI
    // Everyone learns how many global particles exist
    MPI_Bcast(&n_particles_global, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (domain_decomposed) {
        // ========= MpiNaive & BarnesHutMpi: domain-decomposed layout =========
        if (rank == 0) {
            global_particles = sim.particles();  // copy full array (HACC/random/solar)
        }

        counts.resize(processes_count);
        displs.resize(processes_count);

        int base = static_cast<int>(n_particles_global / processes_count);
        int rem  = static_cast<int>(n_particles_global % processes_count);

        for (int r = 0; r < processes_count; ++r) {
            counts[r] = base + (r < rem ? 1 : 0);
        }

        displs[0] = 0;
        for (int r = 1; r < processes_count; ++r) {
            displs[r] = displs[r - 1] + counts[r - 1];
        }

        n_local = counts[rank];
        sim.particles().resize(n_local);

        MPI_Scatterv(
            rank == 0 ? global_particles.data() : nullptr,
            counts.data(), displs.data(), mpi_particle_type,
            sim.particles().data(), n_local, mpi_particle_type,
            0, MPI_COMM_WORLD
        );

        std::printf("Rank %d: local n = %d\n", rank, n_local);
    } else {
        // ========= (fallback replicated layout for future algorithms) =========
        if (rank != 0) {
            sim.particles().resize(n_particles_global);
        }

        MPI_Bcast(
            sim.particles().data(),
            static_cast<int>(n_particles_global),
            mpi_particle_type,
            0, MPI_COMM_WORLD
        );

        std::printf("Rank %d: has full %zu particles\n",
                    rank, n_particles_global);
    }
#else
    // Non-MPI: whole system on single process
    if (n_particles_global == 0)
    {
        std::cerr << "No particles initialized\n";
        return 1;
    }
#endif

    // ----------------------------
    // Performance logger & output
    // ----------------------------
    PerformanceLogger logger(out, RunMode::Run);
#ifdef USE_MPI
    if (rank == 0 && !logger.ok())
        std::cerr << "Warning: could not open output log: " << out << "\n";
#else
    if (!logger.ok())
        std::cerr << "Warning: could not open output log: " << out << "\n";
#endif

    size_t log_step_size = 100;

    // -----------------
    // Checkpoint header
    // -----------------
    if (!disable_output)
    {
#ifdef USE_MPI
    if (rank == 0)
    {
        CheckpointManager* checkpoint = CheckpointManager::getInstance();
        checkpoint->setFilePath("simulation_output.bin");

        SimulationOutputHeader header;
        header.n_particles  = n_particles_global;
        header.target_steps = steps;
        header.passed_steps = 0;
        checkpoint->write_header(header);
        checkpoint->write_masses(global_particles.data(), n_particles_global);
    }
#else
    {
        CheckpointManager* checkpoint = CheckpointManager::getInstance();
        checkpoint->setFilePath("simulation_output.bin");

        SimulationOutputHeader header;
        header.n_particles  = n_particles_global;
        header.target_steps = steps;
        header.passed_steps = 0;
        checkpoint->write_header(header);
        checkpoint->write_masses(sim.particles().data(), n_particles_global);
    }
#endif
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    auto simulation_start = std::chrono::high_resolution_clock::now();

#ifdef USE_MPI
    // Buffer on rank 0 for gathering full system each step
    std::vector<particle_t> global_buffer;
    if (rank == 0 && domain_decomposed)
    {
        global_buffer.resize(n_particles_global);
    }
#endif

    // ---------------
    // Main sim loop
    // ---------------
    for (size_t step = 0; step < steps; ++step)
    {
        sim.step();

#ifdef USE_MPI
        if (domain_decomposed) {
            MPI_Gatherv(
                sim.particles().data(), n_local, mpi_particle_type,
                rank == 0 ? global_buffer.data() : nullptr,
                counts.data(), displs.data(), mpi_particle_type,
                0, MPI_COMM_WORLD
            );

            if (rank == 0)
            {
                if (!disable_output)
                {
                    CheckpointManager* checkpoint = CheckpointManager::getInstance();
                    checkpoint->write_step(global_buffer.data(), n_particles_global);

                    if (step % log_step_size == 0) {
                        std::cout << "Step " << step << " / " << steps << " completed\n";
                    }
                }
            }
        } else {
            // ----- BarnesHutMpi: full data on each rank -----
            if (rank == 0) {
                if (!disable_output)
                {
                    CheckpointManager* checkpoint = CheckpointManager::getInstance();
                    checkpoint->write_step(sim.particles().data(), n_particles_global);

                    if (step % log_step_size == 0) {
                        std::cout << "Step " << step << " / " << steps << " completed\n";
                    }
                }
            }
        }
#else
        if (!disable_output)
        {
            // Sequential
            CheckpointManager* checkpoint = CheckpointManager::getInstance();
            checkpoint->write_step(sim.particles().data(), n_particles_global);

            if (step % log_step_size == 0) {
                std::cout << "Step " << step << " / " << steps << " completed\n";
            }
        }
#endif
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    auto simulation_end = std::chrono::high_resolution_clock::now();
    double simulation_time =
        std::chrono::duration<double>(simulation_end - simulation_start).count();

#ifdef USE_MPI
    if (rank == 0)
    {
        logger.log(n_particles_global, steps, simulation_time, processes_count);
        std::cout << "Simulation completed in " << simulation_time << " seconds\n";
    }
#else
    logger.log(n_particles_global, steps, simulation_time);
    std::cout << "Simulation completed in " << simulation_time << " seconds\n";
#endif

#ifdef USE_MPI
    MPI_Type_free(&mpi_particle_type);
    MPI_Finalize();
#endif

    return 0;
}
