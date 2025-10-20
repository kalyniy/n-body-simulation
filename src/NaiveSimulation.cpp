#include "NaiveSimulation.h"

#include <stdlib.h>
#include <cmath>
#include <vector>

void NaiveSimulation::computeStep(std::vector<particle_t>& particles, const SimParams& params) {
    computeAccelerations_(particles, params);
    integrate_(particles, params);
}

void NaiveSimulation::computeAccelerations_(std::vector<particle_t>& particles, const SimParams& params)
{
    const float eps2 = params.min_r2;
    const float G = params.G;
    const size_t n = particles.size();

    for (auto &p : particles)
        p.acceleration = {0, 0, 0};

    for (size_t i = 0; i < n; ++i)
    {
        particle_t &a = particles[i];

        for (size_t j = i + 1; j < n; ++j)
        {
            particle_t &b = particles[j];
            vector3_t r = b.position - a.position;
            float r2 = r.x * r.x + r.y * r.y + r.z * r.z + eps2;

            float inv_r = 1.0f / std::sqrt(r2);
            float inv_r3 = inv_r * inv_r * inv_r;

            vector3_t acc = r * (G * inv_r3);
            vector3_t a_acc_d = acc * b.mass;
            vector3_t b_acc_d = acc * a.mass;
            a.acceleration += a_acc_d;
            b.acceleration -= b_acc_d;
        }
    }
}

void NaiveSimulation::integrate_(std::vector<particle_t>& particles, const SimParams& params)
{
    const float dt = params.dt;
    for (auto &p : particles)
    {
        p.velocity += p.acceleration * dt;
        p.position += p.velocity * dt;
    }
}