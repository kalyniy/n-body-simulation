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

int WORLD_WIDTH = 300;
int WORLD_HEIGHT = 300;
int WORLD_DEPTH = 300;

float MIN = 0.0001f;

float G = 6.674e-11f;
float dt = 0.01f;

std::vector<particle_t> *g_particles = nullptr;

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
        for (size_t j = i + 1; j < size; j++) // Avoid double calculation and self-interaction
        {
            auto &particle1 = particles->at(i); // Reference to actual particle
            auto &particle2 = particles->at(j); // Reference to actual particle

            auto r = particle2.position - particle1.position;

            // Calculate distance (including z component!)
            auto mag_sq = (r.x * r.x) + (r.y * r.y) + (r.z * r.z);
            auto mag = std::sqrt(mag_sq);

            if (mag > MIN) // Avoid division by zero
            {
                // Gravitational force magnitude: F = G * m1 * m2 / r^2
                float force_mag = G * particle1.mass * particle2.mass / std::max(mag_sq, MIN);

                // Unit vector from particle1 to particle2
                auto unit_r = r * (1.0f / mag);

                // Force vector
                auto force = unit_r * force_mag;

                // Apply Newton's second law: F = ma, so a = F/m
                particle1.acceleration += force * (1.0f / particle1.mass);
                particle2.acceleration += force * (-1.0f / particle2.mass); // Equal and opposite
            }
        }
    }

    // Update velocities and positions using Euler integration
    for (size_t i = 0; i < size; i++)
    {
        auto &particle = particles->at(i);

        // v = v + a * dt
        particle.velocity += particle.acceleration * dt;

        // x = x + v * dt
        particle.position += particle.velocity * dt;

        // Optional: Apply boundary conditions (wrap around or bounce)
        // This prevents particles from flying off to infinity
        if (particle.position.x < 0)
            particle.position.x += WORLD_WIDTH;
        if (particle.position.x > WORLD_WIDTH)
            particle.position.x -= WORLD_WIDTH;
        if (particle.position.y < 0)
            particle.position.y += WORLD_HEIGHT;
        if (particle.position.y > WORLD_HEIGHT)
            particle.position.y -= WORLD_HEIGHT;
        if (particle.position.z < 0)
            particle.position.z += WORLD_DEPTH;
        if (particle.position.z > WORLD_DEPTH)
            particle.position.z -= WORLD_DEPTH;
    }
}

double calculate_elapsed_time(struct timeval start, struct timeval stop)
{
    return (stop.tv_sec + stop.tv_usec * 1e-6) - (start.tv_sec + start.tv_usec * 1e-6);
}

#pragma region OpenGl setup

// Camera controls
float camera_distance = 500.0f;
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
    float mass_normalized = std::min(particle.mass / 1000.0f, 1.0f);
    glColor3f(mass_normalized, 0.5f, 1.0f - mass_normalized);

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

        printf("Generating %lu particles randomly...\n", n);

        generateRandomParticles(g_particles, n);

        printf("Finished generating %lu particles.\n", n);

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
