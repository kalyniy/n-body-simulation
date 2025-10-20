#include "Simulation.h"
#include "renderers/GlutRenderer.h"
#include "NaiveSimulation.h"
#include "BarnesHutSimulation.h"
#include "World.h"

double log_base_8(double x) 
{
    return std::log10(x) / std::log10(8.0);
}

int main(int argc, char **argv)
{
    size_t n_particles = 10000;
   
    SimParams params;
    params.G = 1.0f;
    params.dt = 0.05f; // Smaller timestep for stability!
    params.min_r2 = 2.0f;
    
    auto max_depth = log_base_8((double)n_particles);
    printf("Max depth: %lf\n", max_depth);

    //NBodySimulation sim = NBodySimulation(std::make_unique<NaiveSimulation>(), params);
    NBodySimulation sim = NBodySimulation(std::make_unique<BarnesHutSimulation>(0.5f, 8, max_depth), params);
    
    // Simple choices: default = solar system; pass a folder to load HACC
    if (argc > 1)
    {
        sim.loadHACC(argv[1]);
    }
    else
    {
        //sim.generateClusters(2, n_particles / 2, 300.0f);
        
        sim.generateGalaxyDisk(n_particles, 
            300.0f,   // disk radius
            50.0f     // disk thickness
        );
        
        //sim.generateRandom(n_particles, WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH, 1.0f, 1.0f);
        //sim.setupSolarSystem(WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH);
    }

    GlutRenderer renderer(&argc, argv, 1920, 1080);
    //renderer.attachParticles(&sim.particles());
    renderer.attachSimulation(&sim);

    renderer.onStep([&]()
                    { sim.step(); });

    // Hands off – GLUT runs the loop
    renderer.processEvents();
    return 0;
}