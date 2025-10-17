#include <iostream>
#include <GL/glut.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <algorithm>

#include "Particle.hpp"
#include "DatasetLoader.h"
#include <sys/time.h>

int WORLD_WIDTH = 300 * 2;
int WORLD_HEIGHT = 300 * 2;
int WORLD_DEPTH = 300 * 2;

float MIN = 0.0001f;
float G = 1.0f; // Simulation gravitational constant
float dt = 0.05f;

std::vector<particle_t> *g_particles = nullptr;

void setupSolarSystem(std::vector<particle_t> *particles)
{
    // Clear any existing particles
    particles->clear();
    particles->resize(9);

    // Sun at center - make it much more massive for simulation
    auto &sun = particles->at(0);
    sun.mass = 1000.0f; // Increased mass for better gravitational effect
    sun.position.x = WORLD_WIDTH / 2.0f;
    sun.position.y = WORLD_HEIGHT / 2.0f;
    sun.position.z = WORLD_DEPTH / 2.0f;
    sun.velocity.x = 0.0f;
    sun.velocity.y = 0.0f;
    sun.velocity.z = 0.0f;

    // Mercury
    auto &mercury = particles->at(1);
    mercury.mass = 0.01f; // Scaled masses for simulation
    mercury.position.x = sun.position.x + 60.0f;
    mercury.position.y = sun.position.y;
    mercury.position.z = sun.position.z;
    mercury.velocity.x = 0.0f;
    mercury.velocity.y = 1.2f; // Reduced velocities
    mercury.velocity.z = 0.0f;

    // Venus
    auto &venus = particles->at(2);
    venus.mass = 0.02f;
    venus.position.x = sun.position.x + 90.0f;
    venus.position.y = sun.position.y;
    venus.position.z = sun.position.z;
    venus.velocity.x = 0.0f;
    venus.velocity.y = 1.0f;
    venus.velocity.z = 0.0f;

    // Earth
    auto &earth = particles->at(3);
    earth.mass = 0.02f;
    earth.position.x = sun.position.x + 120.0f;
    earth.position.y = sun.position.y;
    earth.position.z = sun.position.z;
    earth.velocity.x = 0.0f;
    earth.velocity.y = 0.9f;
    earth.velocity.z = 0.0f;

    // Mars
    auto &mars = particles->at(4);
    mars.mass = 0.015f;
    mars.position.x = sun.position.x + 160.0f;
    mars.position.y = sun.position.y;
    mars.position.z = sun.position.z;
    mars.velocity.x = 0.0f;
    mars.velocity.y = 0.75f;
    mars.velocity.z = 0.0f;

    // Jupiter
    auto &jupiter = particles->at(5);
    jupiter.mass = 0.5f; // Much more massive
    jupiter.position.x = sun.position.x + 220.0f;
    jupiter.position.y = sun.position.y;
    jupiter.position.z = sun.position.z;
    jupiter.velocity.x = 0.0f;
    jupiter.velocity.y = 0.6f;
    jupiter.velocity.z = 0.0f;

    // Saturn
    auto &saturn = particles->at(6);
    saturn.mass = 0.3f;
    saturn.position.x = sun.position.x + 280.0f;
    saturn.position.y = sun.position.y;
    saturn.position.z = sun.position.z;
    saturn.velocity.x = 0.0f;
    saturn.velocity.y = 0.5f;
    saturn.velocity.z = 0.0f;

    // Uranus
    auto &uranus = particles->at(7);
    uranus.mass = 0.1f;
    uranus.position.x = sun.position.x;
    uranus.position.y = sun.position.y + 200.0f; // Different axis
    uranus.position.z = sun.position.z;
    uranus.velocity.x = 0.45f; // Velocity in x direction
    uranus.velocity.y = 0.0f;
    uranus.velocity.z = 0.0f;

    // Neptune
    auto &neptune = particles->at(8);
    neptune.mass = 0.1f;
    neptune.position.x = sun.position.x;
    neptune.position.y = sun.position.y - 240.0f; // Opposite side
    neptune.position.z = sun.position.z;
    neptune.velocity.x = -0.4f; // Velocity in negative x direction
    neptune.velocity.y = 0.0f;
    neptune.velocity.z = 0.0f;

    // Add some z-axis variation
    mercury.position.z += 10.0f;
    venus.position.z -= 5.0f;
    mars.position.z += 15.0f;
    jupiter.position.z -= 8.0f;
    saturn.position.z += 12.0f;
    uranus.position.z -= 15.0f;
    neptune.position.z += 8.0f;

    auto set_circular_xy = [&](particle_t &body, const particle_t &sun)
    {
        // Radial vector in the XY plane (ignore z so we keep your small inclination)
        float rx = body.position.x - sun.position.x;
        float ry = body.position.y - sun.position.y;
        float r = std::sqrt(rx * rx + ry * ry);
        if (r < 1e-6f)
            return;

        // Circular speed: v = sqrt(G * M_sun / r)
        float v = std::sqrt(G * sun.mass / r);

        // Tangential unit vector (radial rotated +90° in XY)
        float tx = -ry / r;
        float ty = rx / r;

        body.velocity.x = tx * v;
        body.velocity.y = ty * v;
        body.velocity.z = 0.0f; // keep in-plane; your z position offset gives a tilt
    };

    set_circular_xy(mercury, sun);
    set_circular_xy(venus, sun);
    set_circular_xy(earth, sun);
    set_circular_xy(mars, sun);
    set_circular_xy(jupiter, sun);
    set_circular_xy(saturn, sun);
    set_circular_xy(uranus, sun);
    set_circular_xy(neptune, sun);
}

float randomFloat(int min, int max)
{
    float difference = (float)(max - min);
    float number = min + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (difference)));
    return number;
}

void generateRandomParticles(std::vector<particle_t> *particles, size_t n)
{
    for (size_t i = 0; i < n; i++)
    {
        particle_t particle = {0};

        particle.position = {0};

        particle.position.x = randomFloat(0, WORLD_WIDTH);
        particle.position.y = randomFloat(0, WORLD_HEIGHT);
        particle.position.z = randomFloat(0, WORLD_DEPTH);

        particle.mass = randomFloat(1, 1000);
        particle.acceleration = {0, 0, 0};
        particle.velocity = {0, 0, 0};

        particles->push_back(particle);
    }
}

void computeAccelerations(std::vector<particle_t> *particles)
{
    const float eps2 = 0.0f; // or small softening like 1e-3f if you ever get close passes
    size_t n = particles->size();
    for (size_t i = 0; i < n; ++i)
        particles->at(i).acceleration = {0, 0, 0};

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            auto &a = particles->at(i), &b = particles->at(j);
            auto r = b.position - a.position;
            float r2 = r.x * r.x + r.y * r.y + r.z * r.z + eps2;
            float inv_r = 1.0f / std::sqrt(r2);
            float inv_r3 = inv_r * inv_r * inv_r;
            float s = G * inv_r3;
            // Force per unit mass
            auto acc = r * s;
            a.acceleration += acc * b.mass;
            b.acceleration -= acc * a.mass;
        }
    }
}

void update(std::vector<particle_t> *particles)
{
    auto size = particles->size();

    // Reset accelerations
    for (size_t i = 0; i < size; i++)
    {
        particles->at(i).acceleration = {0, 0, 0};
    }

    // Calculate forces between all particle pairs
    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = i + 1; j < size; j++)
        {
            auto &particle1 = particles->at(i);
            auto &particle2 = particles->at(j);

            auto r = particle2.position - particle1.position;
            auto mag_sq = (r.x * r.x) + (r.y * r.y) + (r.z * r.z);
            auto mag = std::sqrt(mag_sq);

            if (mag > MIN)
            {
                float force_mag = G * particle1.mass * particle2.mass / std::max(mag_sq, MIN);
                auto unit_r = r * (1.0f / mag);
                auto force = unit_r * force_mag;

                particle1.acceleration += force * (1.0f / particle1.mass);
                particle2.acceleration += force * (-1.0f / particle2.mass);
            }
        }
    }

    // Update velocities and positions
    for (size_t i = 0; i < size; i++)
    {
        auto &particle = particles->at(i);

        particle.velocity += particle.acceleration * dt;
        particle.position += particle.velocity * dt;

        // NO BOUNDARY WRAPPING for solar system simulation
        // Let planets orbit freely without artificial boundaries
    }
}

double calculate_elapsed_time(struct timeval start, struct timeval stop)
{
    return (stop.tv_sec + stop.tv_usec * 1e-6) - (start.tv_sec + start.tv_usec * 1e-6);
}

#pragma region OpenGl setup

// Camera controls
float camera_distance = 800.0f;
float camera_angle_x = 0.0f;
float camera_angle_y = 0.0f;
int mouse_x, mouse_y;
bool mouse_dragging = false;

// Animation control
bool animate = true;
int frame_count = 0;

void drawParticle(const particle_t &particle)
{
    glPushMatrix();

    // Translate to particle position
    glTranslatef(particle.position.x - WORLD_WIDTH / 2,
                 particle.position.y - WORLD_HEIGHT / 2,
                 particle.position.z - WORLD_DEPTH / 2);

    // Color based on mass (heavier = redder, lighter = bluer)
    // float mass_normalized = std::min(particle.mass / 1000.0f, 1.0f);
    // glColor3f(mass_normalized, 0.5f, 1.0f - mass_normalized);
    glColor3f(1.0f, 0.5f, 1.0f - 1.0f);

    // Size based on mass
    // float radius = std::max(1.0f, std::sqrt(particle.mass) * 0.1f);
    float radius = 2.0f;

    // Draw sphere
    glutSolidSphere(radius, 8, 6);

    glPopMatrix();
}

void drawBoundingBox()
{
    glColor3f(0.3f, 0.3f, 0.3f);
    glLineWidth(1.0f);

    float half_width = WORLD_WIDTH / 2.0f;
    float half_height = WORLD_HEIGHT / 2.0f;
    float half_depth = WORLD_DEPTH / 2.0f;

    glBegin(GL_LINES);

    // Bottom face
    glVertex3f(-half_width, -half_height, -half_depth);
    glVertex3f(half_width, -half_height, -half_depth);

    glVertex3f(half_width, -half_height, -half_depth);
    glVertex3f(half_width, -half_height, half_depth);

    glVertex3f(half_width, -half_height, half_depth);
    glVertex3f(-half_width, -half_height, half_depth);

    glVertex3f(-half_width, -half_height, half_depth);
    glVertex3f(-half_width, -half_height, -half_depth);

    // Top face
    glVertex3f(-half_width, half_height, -half_depth);
    glVertex3f(half_width, half_height, -half_depth);

    glVertex3f(half_width, half_height, -half_depth);
    glVertex3f(half_width, half_height, half_depth);

    glVertex3f(half_width, half_height, half_depth);
    glVertex3f(-half_width, half_height, half_depth);

    glVertex3f(-half_width, half_height, half_depth);
    glVertex3f(-half_width, half_height, -half_depth);

    // Vertical edges
    glVertex3f(-half_width, -half_height, -half_depth);
    glVertex3f(-half_width, half_height, -half_depth);

    glVertex3f(half_width, -half_height, -half_depth);
    glVertex3f(half_width, half_height, -half_depth);

    glVertex3f(half_width, -half_height, half_depth);
    glVertex3f(half_width, half_height, half_depth);

    glVertex3f(-half_width, -half_height, half_depth);
    glVertex3f(-half_width, half_height, half_depth);

    glEnd();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    // Set up camera
    glTranslatef(0.0f, 0.0f, -camera_distance);
    glRotatef(camera_angle_x, 1.0f, 0.0f, 0.0f);
    glRotatef(camera_angle_y, 0.0f, 1.0f, 0.0f);

    // Draw bounding box
    drawBoundingBox();

    // Draw all particles
    if (g_particles != nullptr)
    {
        for (const auto &particle : *g_particles)
        {
            drawParticle(particle);
        }
    }

    glutSwapBuffers();

    frame_count++;
    if (frame_count % 60 == 0)
    {
        std::cout << "Frame: " << frame_count << std::endl;
    }
}

void reshape(int width, int height)
{
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(45.0, (double)width / (double)height, 1.0, 2000.0);

    glMatrixMode(GL_MODELVIEW);
}

void mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON)
    {
        if (state == GLUT_DOWN)
        {
            mouse_dragging = true;
            mouse_x = x;
            mouse_y = y;
        }
        else
        {
            mouse_dragging = false;
        }
    }
    else if (button == 3) // Mouse wheel up
    {
        camera_distance *= 0.9f;
        glutPostRedisplay();
    }
    else if (button == 4) // Mouse wheel down
    {
        camera_distance *= 1.1f;
        glutPostRedisplay();
    }
}

void motion(int x, int y)
{
    if (mouse_dragging)
    {
        camera_angle_y += (x - mouse_x) * 0.5f;
        camera_angle_x += (y - mouse_y) * 0.5f;

        mouse_x = x;
        mouse_y = y;

        glutPostRedisplay();
    }
}

void keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
    case ' ':
        animate = !animate;
        std::cout << "Animation " << (animate ? "enabled" : "disabled") << std::endl;
        break;
    case 'r':
    case 'R':
        camera_angle_x = 0.0f;
        camera_angle_y = 0.0f;
        camera_distance = 500.0f;
        glutPostRedisplay();
        break;
    case 27: // Escape key
        exit(0);
        break;
    }
}

void idle()
{
    if (animate && g_particles != nullptr)
    {
        update(g_particles);
        glutPostRedisplay();
    }
}

void setupOpenGL(int *argc, char **argv)
{
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("3D Particle Physics Simulation");

    // Set up OpenGL state
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Set up lighting
    GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
    GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat light_ambient[] = {0.3, 0.3, 0.3, 1.0};

    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);

    // Enable color material
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    glClearColor(0.0, 0.0, 0.0, 1.0);

    // Register callbacks
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);

    std::cout << "Controls:" << std::endl;
    std::cout << "  Mouse drag: Rotate camera" << std::endl;
    std::cout << "  Mouse wheel: Zoom in/out" << std::endl;
    std::cout << "  Spacebar: Toggle animation" << std::endl;
    std::cout << "  R: Reset camera" << std::endl;
    std::cout << "  ESC: Exit" << std::endl;
}
#pragma endregion

int main(int argc, char **argv)
{
    srand(static_cast<unsigned>(time(0)));
    // srand(time(NULL));

    g_particles = new std::vector<particle_t>();

    if (argc > 1)
    {
        // auto url = "/home/dan/Desktop/1billionparticles_onesnapshot";
        printf("Trying to load HACC dataset from '%s'\n", argv[1]);
        DatasetLoader::load_hacc_snapshot(g_particles, argv[1]);

        size_t size = g_particles->size();

        std::cout << size << std::endl;
#ifdef DEBUG
        for (size_t i = 0; i < 100; i++)
        {
            particle_t &particle = g_particles->at(i);
            printf("Particle[%ld] x: %f, y: %f, z: %f\n",
                   i,
                   particle.position.x,
                   particle.position.y,
                   particle.position.z);
        }
#endif
    }
    else
    {
        size_t n = 4096;

#ifdef N_PARTICLES
        n = N_PARTICLES;
#endif
        /*
        printf("Generating %lu particles randomly...\n", n);

        generateRandomParticles(g_particles, n);

        printf("Finished generating %lu particles.\n", n);
        */

        setupSolarSystem(g_particles);

        /*
        double avg = 0.0;
        double sum = 0.0;
        double time = 0.0;

        struct timeval start, end;

        for (size_t step = 0; step < g_particles->size(); step++)
        {

            gettimeofday(&start, 0);
            update(g_particles);
            gettimeofday(&end, 0);

            time = calculate_elapsed_time(start, end);

            sum += time;

            avg = sum / (step + 1);

            printf("Avg step time: %lf\n", avg);
#ifdef DEBUG
            if (step % 100 == 0)
            {
                auto &particle = g_particles->at(0);
                printf("Particle[%uld] x: %f, y: %f, z: %f\n",
                       0,
                       particle.position.x,
                       particle.position.y,
                       particle.position.z);
            }
#endif
        }
        */
    }

    setupOpenGL(&argc, argv);
    glutMainLoop();

    delete g_particles;
    return 0;
}
