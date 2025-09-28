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

int WIDTH = 300;
int LENGTH = 300;
int HEIGHT = 300;

float MIN = 0.0001f;

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
        std::cout << i << std::endl;
        particle_t particle = {0};

        particle.position = {0};

        particle.position.x = randomFloat(0, WIDTH);
        particle.position.y = randomFloat(0, LENGTH);
        particle.position.z = randomFloat(0, HEIGHT);

        particle.mass = randomFloat(0, INT32_MAX);
        particle.acceleration = 0;

        particles->push_back(particle);
    }
}

void update(std::vector<particle_t> *particles)
{
    auto size = particles->size();
    for (size_t i = 0; i < size; i++)
    {
        auto particle = particles->at(i);
        auto p1 = particle.position;

        for (size_t j = 0; j < size; j++)
        {
            if (i == j)
            {
                continue;
            }

            auto particle2 = particles->at(j);
            auto p2 = particle2.position;
            auto m2 = particle2.mass;

            auto r = p2 - p1;

            auto mag_sq = (r.x * r.x) + (r.y * r.y);
            auto mag = std::sqrt(mag_sq);

            float part1 = (m2 / (std::max(mag_sq, MIN) * mag));

            auto a1 = r * part1;

            particle.acceleration = 1;
        }
    }
}

int main(int argc, char **argv)
{
    srand(static_cast<unsigned>(time(0)));
    // srand(time(NULL));

    std::vector<particle_t> *particles = new std::vector<particle_t>();

    if (argc > 1)
    {
        // auto url = "/home/dan/Desktop/1billionparticles_onesnapshot";

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

        generateRandomParticles(particles, n);

        std::cout << particles->size() << std::endl;

        for (size_t i = 0; i < particles->size(); i++)
        {
            sleep(1);
            update(particles);
            std::cout << "Updated particles" << std::endl;
        }
    }

    setupOpenGl(&argc, argv);

    delete particles;
    return 0;
}
