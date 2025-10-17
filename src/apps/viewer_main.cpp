#include "Simulation.h"
#include "renderers/GlutRenderer.h"

int main(int argc, char **argv)
{
    SimParams params;
    params.G = 1.0f;
    params.dt = 0.05f;
    NBodySimulation sim(params);

    // Simple choices: default = solar system; pass a folder to load HACC
    if (argc > 1)
        sim.loadHACC(argv[1]);
    else
        sim.setupSolarSystem(600, 600, 600);

    GlutRenderer renderer(&argc, argv, 1280, 800);
    renderer.attachParticles(&sim.particles());
    renderer.onStep([&]()
                    { sim.step(); });

    // Hands off – GLUT runs the loop
    renderer.processEvents();
    return 0;
}
