#include "Simulation.h"
#include "renderers/GlutRenderer.h"
#include "NaiveSimulation.h"
#include "World.h"
#include "CheckpointManager.h"
#include <thread>
#include <chrono>
#include <fstream>

struct PlaybackController
{
    size_t current_step = 0;
    size_t total_steps = 0;
    bool playing = true;
    bool loop = false;
    int speed = 1;
    int frame_delay_ms = 16; // ~60 FPS
    CheckpointManager* checkpoint;
    NBodySimulation* sim;
    float* positions;
    size_t n_particles;
    
    void loadStep(size_t step)
    {
        if (step >= total_steps) return;
        
        size_t n_read = checkpoint->read_step(positions, step);
        
        if (n_read > 0)
        {
            // Update particle positions directly
            for (size_t i = 0; i < n_read; ++i)
            {
                sim->particles()[i].position.x = positions[i * 3 + 0];
                sim->particles()[i].position.y = positions[i * 3 + 1];
                sim->particles()[i].position.z = positions[i * 3 + 2];
            }
            current_step = step;
            
            // Force update of render buffer
            sim->updateRenderBuffer();
        }
    }
    
    void advanceStep()
    {
        if (!playing) return;
        
        for (int s = 0; s < speed; ++s)
        {
            current_step++;
            
            if (current_step >= total_steps)
            {
                if (loop)
                {
                    current_step = 0;
                    std::cout << "Looping back to start" << std::endl;
                }
                else
                {
                    current_step = total_steps - 1;
                    playing = false;
                    std::cout << "Playback complete!" << std::endl;
                    return;
                }
            }
        }
        
        loadStep(current_step);
        
        if (current_step % 100 == 0)
        {
            std::cout << "Step: " << current_step << " / " << total_steps << std::endl;
        }
    }
    
    void jumpToStep(size_t step)
    {
        if (step < total_steps)
        {
            loadStep(step);
            std::cout << "Jumped to step: " << current_step << std::endl;
        }
    }
    
    void reset()
    {
        loadStep(0);
        playing = false;
        std::cout << "Reset to beginning (paused)" << std::endl;
    }
    
    void increaseSpeed()
    {
        if (speed < 10)
        {
            speed++;
            std::cout << "Speed: " << speed << "x" << std::endl;
        }
    }
    
    void decreaseSpeed()
    {
        if (speed > 1)
        {
            speed--;
            std::cout << "Speed: " << speed << "x" << std::endl;
        }
    }
    
    void toggleLoop()
    {
        loop = !loop;
        std::cout << "Loop: " << (loop ? "ON" : "OFF") << std::endl;
    }
    
    void togglePlayPause()
    {
        playing = !playing;
        std::cout << (playing ? "Playing" : "Paused") 
                  << " - Step: " << current_step << " / " << total_steps << std::endl;
    }
    
    void skipForward()
    {
        size_t skip = total_steps / 10; // Skip 10%
        if (skip < 10) skip = 10;
        jumpToStep(current_step + skip);
    }
    
    void skipBackward()
    {
        size_t skip = total_steps / 10;
        if (skip < 10) skip = 10;
        if (current_step > skip)
            jumpToStep(current_step - skip);
        else
            jumpToStep(0);
    }
};

int main(int argc, char **argv)
{
    std::string checkpoint_file = "simulation_output.bin";
    int playback_speed = 1;
    int frame_delay = 16; // milliseconds between frames
    
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (a == "--file" && i + 1 < argc)
            checkpoint_file = argv[++i];
        else if (a == "--speed" && i + 1 < argc)
            playback_speed = std::stoi(argv[++i]);
        else if (a == "--fps" && i + 1 < argc)
            frame_delay = 1000 / std::stoi(argv[++i]);
    }
    
    // Check if file exists
    std::ifstream test(checkpoint_file);
    if (!test.good())
    {
        std::cerr << "Error: Checkpoint file not found: " << checkpoint_file << std::endl;
        std::cerr << "Please run the MPI simulation first to generate the data." << std::endl;
        return 1;
    }
    test.close();
    
    // Initialize checkpoint manager
    CheckpointManager* checkpoint = CheckpointManager::getInstance();
    checkpoint->setFilePath(checkpoint_file);
    
    // Read header
    SimulationOutputHeader header = checkpoint->read_header();
    
    // Verify simulation is complete
    if (header.passed_steps < header.target_steps)
    {
        std::cout << "Warning: Simulation appears incomplete.\n"
                  << "  Expected steps: " << header.target_steps << "\n"
                  << "  Available steps: " << header.passed_steps << "\n"
                  << "  Proceeding with available data...\n" << std::endl;
    }
    
    std::cout << "Loaded simulation data:\n"
              << "  Particles: " << header.n_particles << "\n"
              << "  Total steps: " << header.passed_steps << "\n"
              << "  File: " << checkpoint_file << std::endl;
    
    // Create simulation object (just for holding particle data)
    SimParams params;
    NBodySimulation sim = NBodySimulation(std::make_unique<NaiveSimulation>(), params);
    sim.particles().resize(header.n_particles);
    
    // Initialize particles with default values
    for (size_t i = 0; i < header.n_particles; ++i)
    {
        sim.particles()[i].mass = 1.0f;
        sim.particles()[i].velocity = {0, 0, 0};
        sim.particles()[i].acceleration = {0, 0, 0};
    }
    
    // Allocate position buffer
    float* positions = new float[header.n_particles * 3];
    
    // Create playback controller
    PlaybackController* controller = new PlaybackController{
        0,
        header.passed_steps,
        true,  // Start playing
        false, // Loop off
        playback_speed,
        frame_delay,
        checkpoint,
        &sim,
        positions,
        header.n_particles
    };
    
    // Load first step
    std::cout << "Loading initial step..." << std::endl;
    controller->loadStep(0);
    std::cout << "First particle position: (" 
              << sim.particles()[0].position.x << ", "
              << sim.particles()[0].position.y << ", "
              << sim.particles()[0].position.z << ")" << std::endl;
    
    // Initialize renderer
    GlutRenderer renderer(&argc, argv, 1920, 1080);
    renderer.attachSimulation(&sim);
    
    // Set up playback step callback
    renderer.onStep([controller]()
    {
        controller->advanceStep();
        // Optional: add frame delay to control playback speed
        std::this_thread::sleep_for(std::chrono::milliseconds(controller->frame_delay_ms));
    });
    
    // Set up keyboard controls
    renderer.onKeyboard([controller](unsigned char key, int x, int y)
    {
        switch (key)
        {
        case 'l':
        case 'L':
            controller->toggleLoop();
            break;
        case '+':
        case '=':
            controller->increaseSpeed();
            break;
        case '-':
        case '_':
            controller->decreaseSpeed();
            break;
        case '0':
            controller->reset();
            break;
        case 'p':
        case 'P':
            controller->togglePlayPause();
            break;
        case '[':
            controller->skipBackward();
            break;
        case ']':
            controller->skipForward();
            break;
        }
    });
    
    std::cout << "\n=== Playback Controls ===\n"
              << "  Space/P: Play/Pause\n"
              << "  L: Toggle loop\n"
              << "  +/-: Increase/Decrease speed\n"
              << "  [/]: Skip backward/forward 10%\n"
              << "  0: Reset to beginning\n"
              << "  R: Reset camera\n"
              << "  Mouse drag: Rotate camera\n"
              << "  Mouse wheel: Zoom\n"
              << "  ESC: Quit\n" 
              << "========================\n" << std::endl;
    
    // Run the renderer
    renderer.processEvents();
    
    // Cleanup
    delete[] positions;
    delete controller;
    
    return 0;
}