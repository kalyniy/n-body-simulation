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

void displayMe(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POLYGON);
    glVertex3f(0.5, 0.0, 0.5);
    glVertex3f(0.5, 0.0, 0.0);
    glVertex3f(0.0, 0.5, 0.0);
    glVertex3f(0.0, 0.0, 0.5);
    glEnd();
    glFlush();
}

void setupOpenGl(int *argc, char **argv)
{
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(400, 300);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Hello world!");
    glutDisplayFunc(displayMe);
    glutMainLoop();
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

        particle.mass = randomFloat(0, INT32_MAX);
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

int main(int argc, char **argv)
{
    srand(static_cast<unsigned>(time(0)));
    // srand(time(NULL));

    std::vector<particle_t> *particles = new std::vector<particle_t>();

    if (argc > 1)
    {
        // auto url = "/home/dan/Desktop/1billionparticles_onesnapshot";
        printf("Trying to load HACC dataset from '%s'\n", argv[1]);
        DatasetLoader::load_hacc_snapshot(particles, argv[1]);

        size_t size = particles->size();

        std::cout << size << std::endl;
#ifdef DEBUG
        for (size_t i = 0; i < 100; i++)
        {
            particle_t &particle = particles->at(i);
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

        generateRandomParticles(particles, n);

        std::cout << particles->size() << std::endl;

        double avg = 0.0;
        double sum = 0.0;
        double time = 0.0;

        struct timeval start, end;

        for (size_t step = 0; step < particles->size(); step++)
        {

            gettimeofday(&start, 0);
            update(particles);
            gettimeofday(&end, 0);

            time = calculate_elapsed_time(start, end);

            sum += time;

            avg = sum / (step + 1);

            printf("Avg step time: %lf\n", avg);
#ifdef DEBUG
            if (step % 100 == 0)
            {
                auto &particle = particles->at(0);
                printf("Particle[%uld] x: %f, y: %f, z: %f\n",
                       0,
                       particle.position.x,
                       particle.position.y,
                       particle.position.z);
            }
#endif
        }
    }

    setupOpenGl(&argc, argv);

    delete particles;
    return 0;
}
