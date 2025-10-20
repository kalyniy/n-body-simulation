#include "Simulation.h"
#include "renderers/GlutRenderer.h"
#include "NaiveSimulation.h"
#include "BarnesHutSimulation.h"
#include "World.h"

int main(int argc, char **argv)
{
    SimParams params;
    params.G = 1.0f;
    params.dt = 0.05f;
    
    //NBodySimulation sim = NBodySimulation(std::make_unique<NaiveSimulation>(), params);
    NBodySimulation sim = NBodySimulation(std::make_unique<BarnesHutSimulation>(0.5f), params);
    
    // Simple choices: default = solar system; pass a folder to load HACC
    if (argc > 1)
    {
        sim.loadHACC(argv[1]);
    }
    else
    {
        sim.generateRandom(100000, WORLD_WIDTH, WORLD_HEIGHT, WORLD_DEPTH, 1.0f, 1.0f);
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