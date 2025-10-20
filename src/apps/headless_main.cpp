#include <chrono>
#include <iostream>
#include "Simulation.h"
#include "PerformanceLogger.hpp"
#include "NaiveSimulation.h"
#include "BarnesHutSimulation.h"
#include "World.h"

int main(int argc, char **argv)
{
    SimParams params;
    params.G = 1.0f;
    params.dt = 0.05f;

    NBodySimulation sim = NBodySimulation(std::make_unique<BarnesHutSimulation>(0.8f), params);

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
        std::cout << "Generating " << randomN << " random particles\n";
        sim.generateRandom(randomN, WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);
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
    size_t sizes[] = {100, 1000, 2000, 5000, 10000};
    size_t sizes_count = sizeof(sizes) / sizeof(sizes[0]);

    //for (size_t i = 0; i < sizes_count; i++)
    //{
        const size_t n_particles = randomN; //sizes[i];
        std::cout << "Generating " << n_particles << " random particles\n";

        //sim.generateRandom(n_particles, WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);

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
    //}

    return 0;
}