#pragma once
#ifdef __APPLE__
#define MACOS
#endif

#ifdef MACOS
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <vector>
#include <functional>
#include "Renderer.h"
//#include "Particle.hpp"
#include "Simulation.h"
// Minimal GLUT renderer that can "drive" a simulation via a callback.
class GlutRenderer : public Renderer
{
public:
    using StepFn = std::function<void()>; // called each idle frame to advance sim
    using KeyboardFn = std::function<void(unsigned char, int, int)>;
    void onKeyboard(KeyboardFn fn) { keyboard_fn_ = std::move(fn); }

    GlutRenderer(int *pargc, char **argv, int w, int h);
    ~GlutRenderer() override = default;

    //void attachParticles(const std::vector<particle_t> *particles);
    void attachSimulation(NBodySimulation* sim);
    void onStep(StepFn fn); // optional: sim->step() provided from app

    void draw(const std::vector<particle_t> &particles) override; // not used directly
    void processEvents() override;                                // GLUT takes over
    bool shouldClose() const override { return should_close_; }

private:
    static GlutRenderer *instance_;
    static void sDisplay();
    static void sReshape(int w, int h);
    static void sMouse(int button, int state, int x, int y);
    static void sMotion(int x, int y);
    static void sKeyboard(unsigned char key, int x, int y);
    static void sIdle();

    void display_();
    void reshape_(int w, int h);
    void mouse_(int button, int state, int x, int y);
    void motion_(int x, int y);
    void keyboard_(unsigned char key, int x, int y);
    void idle_();

private:
    NBodySimulation* simulation_ = nullptr; 
    //const std::vector<particle_t> *particles_ = nullptr;
    StepFn step_fn_{};
    KeyboardFn keyboard_fn_{};

    float camera_distance_ = 800.f;
    float camera_angle_x_ = 0.f;
    float camera_angle_y_ = 0.f;
    int mouse_x_ = 0, mouse_y_ = 0;
    bool mouse_dragging_ = false;
    bool animate_ = true;
    bool should_close_ = false;
};
