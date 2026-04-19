#include "NaiveSimulation.h"
#include <cmath>

void NaiveSimulation::computeStep(ParticleSoA& p, const SimParams& params) {
    computeAccelerations_(p, params);
    integrate_(p, params);
}

void NaiveSimulation::computeAccelerations_(ParticleSoA& p, const SimParams& params)
{
    const float eps2 = params.min_r2;
    const float G = params.G;
    const size_t n = p.size();

    p.zeroAccelerations();

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            float rx = p.pos_x[j] - p.pos_x[i];
            float ry = p.pos_y[j] - p.pos_y[i];
            float rz = p.pos_z[j] - p.pos_z[i];
            float r2 = rx * rx + ry * ry + rz * rz + eps2;

            float inv_r = 1.0f / std::sqrt(r2);
            float inv_r3 = inv_r * inv_r * inv_r;

            float ax = G * inv_r3 * rx;
            float ay = G * inv_r3 * ry;
            float az = G * inv_r3 * rz;

            p.acc_x[i] += ax * p.mass[j];
            p.acc_y[i] += ay * p.mass[j];
            p.acc_z[i] += az * p.mass[j];

            p.acc_x[j] -= ax * p.mass[i];
            p.acc_y[j] -= ay * p.mass[i];
            p.acc_z[j] -= az * p.mass[i];
        }
    }
}

void NaiveSimulation::integrate_(ParticleSoA& p, const SimParams& params)
{
    const float dt = params.dt;
    for (size_t i = 0; i < p.size(); ++i)
    {
        p.vel_x[i] += p.acc_x[i] * dt;
        p.vel_y[i] += p.acc_y[i] * dt;
        p.vel_z[i] += p.acc_z[i] * dt;

        p.pos_x[i] += p.vel_x[i] * dt;
        p.pos_y[i] += p.vel_y[i] * dt;
        p.pos_z[i] += p.vel_z[i] * dt;
    }
}