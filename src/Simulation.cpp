#include "Simulation.h"
#include "DatasetLoader.h"
#include <random>
#include <algorithm>
#include <cmath>

void NBodySimulation::clear() { particles_.clear(); }
void NBodySimulation::reserve(std::size_t n) { particles_.reserve(n); }

void NBodySimulation::loadHACC(const std::string &dir)
{
    particles_.clear();
    DatasetLoader::load_hacc_snapshot(&particles_, dir);
}

void NBodySimulation::generateRandom(std::size_t n, int W, int H, int D, float minM, float maxM)
{
    particles_.clear();
    particles_.reserve(n);
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> px(0, (float)W), py(0, (float)H), pz(0, (float)D);
    std::uniform_real_distribution<float> pm(minM, maxM);

    for (size_t i = 0; i < n; ++i)
    {
        particle_t p{};
        p.position = {px(rng), py(rng), pz(rng)};
        p.velocity = {0, 0, 0};
        p.acceleration = {0, 0, 0};
        p.mass = pm(rng);
        particles_.push_back(p);
    }
}

void NBodySimulation::setupSolarSystem(int W, int H, int D)
{
    particles_.clear();
    particles_.resize(9);

    auto &sun = particles_[0];
    sun.mass = 1000.0f;
    sun.position = {(float)W / 2.f, (float)H / 2.f, (float)D / 2.f};
    sun.velocity = {0, 0, 0};

    auto set_body = [&](int idx, float rxy, float zoff)
    {
        auto &b = particles_[idx];
        b.mass = 0.02f;
        b.position = {sun.position.x + rxy, sun.position.y, sun.position.z + zoff};
        b.velocity = {0, 0, 0};
    };

    set_body(1, 60.f, +10.f); // Mercury-ish
    set_body(2, 90.f, -5.f);
    set_body(3, 120.f, 0.f); // Earth-ish
    set_body(4, 160.f, +15.f);
    set_body(5, 220.f, -8.f);
    particles_[5].mass = 0.5f; // Jupiter-ish
    set_body(6, 280.f, +12.f);
    particles_[6].mass = 0.3f;
    set_body(7, 0.f, -15.f);
    particles_[7].position = {sun.position.x, sun.position.y + 200.f, sun.position.z - 15.f};
    particles_[7].mass = 0.1f;
    set_body(8, 0.f, +8.f);
    particles_[8].position = {sun.position.x, sun.position.y - 240.f, sun.position.z + 8.f};
    particles_[8].mass = 0.1f;

    auto set_circular_xy = [&](particle_t &body, const particle_t &s)
    {
        float rx = body.position.x - s.position.x;
        float ry = body.position.y - s.position.y;
        float r = std::sqrt(rx * rx + ry * ry);
        if (r < 1e-6f)
            return;
        float v = std::sqrt(params_.G * s.mass / r);
        float tx = -ry / r, ty = rx / r; // tangential
        body.velocity = {tx * v, ty * v, 0.0f};
    };

    for (int i = 1; i < 9; ++i)
        set_circular_xy(particles_[i], sun);
}

/**
 * @brief Resets
 *
 */
void reset()
{
}

/**
 * @brief Computes Accelerations in O(n^2).
 * @brief 20 floating operations
 */
void NBodySimulation::computeAccelerations_()
{
    const float eps2 = params_.min_r2;
    const float G = params_.G;
    const size_t n = particles_.size();

    for (auto &p : particles_)
        p.acceleration = {0, 0, 0};

    for (size_t i = 0; i < n; ++i)
    {
        particle_t &a = particles_[i];
        auto massA = a.mass;

        // vector3_t acc_i = {0, 0, 0}; // Accumulator for particle i
        for (size_t j = i + 1; j < n; ++j)
        {
            particle_t &b = particles_[j];
            auto massB = b.mass;
            vector3_t r = b.position - a.position;
            float r2 = r.x * r.x + r.y * r.y + r.z * r.z + eps2;

            float inv_r = 1.0f / sqrt(r2); // std::sqrt(r2);
            float inv_r3 = inv_r * inv_r * inv_r;

            vector3_t acc = r * (G * inv_r3);
            auto a_acc_d = acc * b.mass;
            auto b_acc_d = acc * a.mass;
            a.acceleration += a_acc_d;
            b.acceleration -= b_acc_d;
        }
    }
}

void NBodySimulation::integrate_()
{
    const float dt = params_.dt;
    for (auto &p : particles_)
    {
        p.velocity += p.acceleration * dt;
        p.position += p.velocity * dt;
    }
}

void NBodySimulation::step()
{
    computeAccelerations_();
    integrate_();
}