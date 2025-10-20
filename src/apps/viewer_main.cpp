#include "Simulation.h"
#include "renderers/GlutRenderer.h"
#include "NaiveSimulation.h"
#include "BarnesHutSimulation.h"

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
        sim.generateRandom(10000, 600, 600, 600, 1.0f, 1.0f);
        //sim.setupSolarSystem(600, 600, 600);
    }

    GlutRenderer renderer(&argc, argv, 1280, 800);
    renderer.attachParticles(&sim.particles());
    renderer.onStep([&]()
                    { sim.step(); });

    // Hands off – GLUT runs the loop
    renderer.processEvents();
    return 0;
}