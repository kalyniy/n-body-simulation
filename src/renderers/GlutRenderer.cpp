#include "renderers/GlutRenderer.h"
#include "Simulation.h"
#include "Utils.hpp"
#include <iostream>
#include <algorithm>
#include <time.h>
#include <sys/time.h>

GlutRenderer *GlutRenderer::instance_ = nullptr;

GlutRenderer::GlutRenderer(int *pargc, char **argv, int w, int h)
{
    instance_ = this;
    glutInit(pargc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(w, h);
    glutCreateWindow("N-Body Viewer (GLUT)");

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    GLfloat light_position[] = {1.f, 1.f, 1.f, 0.f};
    GLfloat light_diffuse[] = {1.f, 1.f, 1.f, 1.f};
    GLfloat light_ambient[] = {0.3f, 0.3f, 0.3f, 1.f};
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glClearColor(0, 0, 0, 1);

    glutDisplayFunc(GlutRenderer::sDisplay);
    glutReshapeFunc(GlutRenderer::sReshape);
    glutMouseFunc(GlutRenderer::sMouse);
    glutMotionFunc(GlutRenderer::sMotion);
    glutKeyboardFunc(GlutRenderer::sKeyboard);
    glutIdleFunc(GlutRenderer::sIdle);

    std::cout << "Controls:\n"
              << "  Mouse drag: rotate\n"
              << "  Wheel: zoom\n"
              << "  Space: toggle animation\n"
              << "  R: reset camera\n"
              << "  ESC: quit\n";
}

//void GlutRenderer::attachParticles(const std::vector<particle_t> *p) { particles_ = p; }

void GlutRenderer::attachSimulation(NBodySimulation* sim)
{
    simulation_ = sim;
}

void GlutRenderer::onStep(StepFn fn) { step_fn_ = std::move(fn); }

void GlutRenderer::processEvents()
{
    // GLUT enters its own loop; this call never returns.
    glutMainLoop();
}

void GlutRenderer::draw(const std::vector<particle_t> &)
{
    // Not used; GLUT uses callbacks instead.
}

void GlutRenderer::sDisplay()
{
    if (instance_)
        instance_->display_();
}
void GlutRenderer::sReshape(int w, int h)
{
    if (instance_)
        instance_->reshape_(w, h);
}
void GlutRenderer::sMouse(int b, int s, int x, int y)
{
    if (instance_)
        instance_->mouse_(b, s, x, y);
}
void GlutRenderer::sMotion(int x, int y)
{
    if (instance_)
        instance_->motion_(x, y);
}
void GlutRenderer::sKeyboard(unsigned char k, int x, int y)
{
    if (instance_)
        instance_->keyboard_(k, x, y);
}
void GlutRenderer::sIdle()
{
    if (instance_)
        instance_->idle_();
}

void GlutRenderer::display_()
{
    if (!simulation_) return;
    
    const auto& render_particles = simulation_->getRenderBuffer();
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    
    glTranslatef(0.f, 0.f, -camera_distance_);
    glRotatef(camera_angle_x_, 1.f, 0.f, 0.f);
    glRotatef(camera_angle_y_, 0.f, 1.f, 0.f);
    
    // Debug output
    static int frame_count = 0;
    if (frame_count % 60 == 0 && !render_particles.empty())
    {
        std::cout << "Rendering " << render_particles.size() << " particles" << std::endl;
        std::cout << "First particle: (" 
                  << render_particles[0].position.x << ", "
                  << render_particles[0].position.y << ", "
                  << render_particles[0].position.z << ")" << std::endl;
    }
    frame_count++;
    
    // Draw particles
    if (!render_particles.empty())
    {
        for (const auto &p : render_particles)
        {
            glPushMatrix();
            glTranslatef(p.position.x, p.position.y, p.position.z);
            glColor3f(1.f, 0.5f, 0.f);
            glutSolidSphere(2.f, 8, 6);
            glPopMatrix();
        }
    }
    else
    {
        std::cout << "Warning: No particles to render!" << std::endl;
    }
    
    glutSwapBuffers();
}

void GlutRenderer::reshape_(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (double)w / (double)h, 1.0, 5000.0);
    glMatrixMode(GL_MODELVIEW);
}

void GlutRenderer::mouse_(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON)
    {
        mouse_dragging_ = (state == GLUT_DOWN);
        mouse_x_ = x;
        mouse_y_ = y;
    }
    else if (button == 3 && state == GLUT_DOWN)
    { // wheel up
        camera_distance_ *= 0.9f;
        glutPostRedisplay();
    }
    else if (button == 4 && state == GLUT_DOWN)
    { // wheel down
        camera_distance_ *= 1.1f;
        glutPostRedisplay();
    }
}

void GlutRenderer::motion_(int x, int y)
{
    if (!mouse_dragging_)
        return;
    camera_angle_y_ += (x - mouse_x_) * 0.5f;
    camera_angle_x_ += (y - mouse_y_) * 0.5f;
    mouse_x_ = x;
    mouse_y_ = y;
    glutPostRedisplay();
}
/*
void GlutRenderer::keyboard_(unsigned char key, int, int)
{
    switch (key)
    {
    case ' ':
        animate_ = !animate_;
        break;
    case 'r':
    case 'R':
        camera_angle_x_ = camera_angle_y_ = 0.f;
        camera_distance_ = 800.f;
        break;
    case 27:
        should_close_ = true;
        std::exit(0);
        break;
    }
}
*/

void GlutRenderer::idle_()
{
    if (animate_ && step_fn_)
        step_fn_();
    glutPostRedisplay();
}

void GlutRenderer::keyboard_(unsigned char key, int x, int y)
{
    // Call custom keyboard handler first if set
    if (keyboard_fn_)
    {
        keyboard_fn_(key, x, y);
    }
    
    // Default controls
    switch (key)
    {
    case ' ':
        animate_ = !animate_;
        std::cout << (animate_ ? "Playing" : "Paused") << std::endl;
        break;
    case 'r':
    case 'R':
        camera_angle_x_ = camera_angle_y_ = 0.f;
        camera_distance_ = 800.f;
        std::cout << "Camera reset" << std::endl;
        break;
    case 27: // ESC
        should_close_ = true;
        std::exit(0);
        break;
    }
}