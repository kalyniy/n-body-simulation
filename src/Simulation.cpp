#include "Simulation.h"
#include "DatasetLoader.h"
#include "Random.hpp"
#include <random>
#include <algorithm>
#include <cmath>
#include <iostream>

void NBodySimulation::clear() { particles_.clear(); }
void NBodySimulation::reserve(std::size_t n) { particles_.reserve(n); }

void NBodySimulation::loadHACC(const std::string &dir)
{
    particles_.clear();
    DatasetLoader::load_hacc_snapshot(&particles_, dir);
}

void NBodySimulation::generateRandom(std::size_t n, int W, int H, int D, float minM, float maxM)
{
    std::cout << "Generating " << n << " random particles\n";

    particles_.clear();
    particles_.reserve(n);
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> px(0, (float)W), py(0, (float)H), pz(0, (float)D);
    std::uniform_real_distribution<float> pm(minM, maxM);

    for (size_t i = 0; i < n; ++i)
        particles_.push_back(px(rng), py(rng), pz(rng), 0, 0, 0, pm(rng));
}

void NBodySimulation::setupSolarSystem(int W, int H, int D)
{
    particles_.clear();
    particles_.resize(9);

    float cx = W/2.f, cy = H/2.f, cz = D/2.f;
    // Sun
    particles_.pos_x[0] = cx;
    particles_.pos_y[0] = cy;
    particles_.pos_z[0] = cz;
    particles_.vel_x[0] = 0;
    particles_.vel_y[0] = 0;
    particles_.vel_z[0] = 0;
    particles_.mass[0] = 1000.f;

    auto set_body = [&](int idx, float rxy, float zoff, float m) {
        particles_.pos_x[idx] = cx+rxy;
        particles_.pos_y[idx] = cy;
        particles_.pos_z[idx] = cz+zoff;
        particles_.vel_x[idx] = 0;
        particles_.vel_y[idx] = 0;
        particles_.vel_z[idx] = 0;
        particles_.mass[idx] = m;
    };

    set_body(1, 60.f, +10.f, 0.02f);
    set_body(2, 90.f, -5.f, 0.02f);
    set_body(3, 120.f, 0.f, 0.02f);
    set_body(4, 160.f, +15.f, 0.02f);
    set_body(5, 220.f, -8.f, 0.5f);
    set_body(6, 280.f, +12.f, 0.3f);
    // Two bodies at different y offsets
    particles_.pos_x[7]=cx; particles_.pos_y[7]=cy+200.f; particles_.pos_z[7]=cz-15.f;
    particles_.vel_x[7]=0; particles_.vel_y[7]=0; particles_.vel_z[7]=0; particles_.mass[7]=0.1f;
    particles_.pos_x[8]=cx; particles_.pos_y[8]=cy-240.f; particles_.pos_z[8]=cz+8.f;
    particles_.vel_x[8]=0; particles_.vel_y[8]=0; particles_.vel_z[8]=0; particles_.mass[8]=0.1f;

    // Set circular orbits
    for (int i = 1; i < 9; ++i) {
        float rx = particles_.pos_x[i] - particles_.pos_x[0];
        float ry = particles_.pos_y[i] - particles_.pos_y[0];
        float r = std::sqrt(rx*rx + ry*ry);
        if (r < 1e-6f) continue;
        float v = std::sqrt(params_.G * particles_.mass[0] / r);
        float tx = -ry/r, ty = rx/r;
        particles_.vel_x[i] = tx*v;
        particles_.vel_y[i] = ty*v;
        particles_.vel_z[i] = 0.f;
    }
}

void NBodySimulation::generateGalaxyDisk(int n_particles, float radius, float thickness)
{
    particles_.clear();
    particles_.reserve(n_particles);

    const float central_mass = 1000000.0f;
    particles_.push_back(0, 0, 0, 0, 0, 0, central_mass);

    for (int i = 1; i < n_particles; ++i) {
        float r = radius * std::sqrt(uniformRandom(0.0f, 1.0f));
        float theta = uniformRandom(0.0f, 2.0f * M_PI);
        float z = gaussianRandom(0.0f, thickness * 0.3f);

        float px = r * std::cos(theta);
        float pz = r * std::sin(theta);

        float v_orbital = std::sqrt(params_.G * central_mass / r);
        float v_d = 0.1f * v_orbital;

        float vx = -v_orbital * std::sin(theta) + gaussianRandom(0.f, v_d);
        float vy = gaussianRandom(0.f, v_d * 0.5f);
        float vz =  v_orbital * std::cos(theta) + gaussianRandom(0.f, v_d);

        particles_.push_back(px, z, pz, vx, vy, vz, uniformRandom(0.5f, 2.0f));
    }

    float disk_volume = M_PI * radius * radius * thickness;
    float avg_spacing = std::cbrt(disk_volume / n_particles);
    float softening = avg_spacing * 0.1f;
    params_.min_r2 = softening * softening;

    std::cout << "Disk: avg_spacing = " << avg_spacing
              << ", softening = " << softening
              << ", min_r2 = " << params_.min_r2 << "\n";
}

void NBodySimulation::generateClusters(int n_clusters, int particles_per_cluster, float cluster_separation)
{
    particles_.clear();
    particles_.reserve(n_clusters * particles_per_cluster);

    for (int c = 0; c < n_clusters; ++c) {
        float angle = (2.0f * M_PI * c) / n_clusters;
        float ccx = cluster_separation * std::cos(angle);
        float ccy = uniformRandom(-50.0f, 50.0f);
        float ccz = cluster_separation * std::sin(angle);

        float speed = 5.0f;
        float cvx = -speed * std::cos(angle);
        float cvz = -speed * std::sin(angle);

        float cluster_radius = 100.0f;
        for (int i = 0; i < particles_per_cluster; ++i) {
            auto gv = gaussianRandomVec3(0.0f, cluster_radius);
            auto rv = gaussianRandomVec3(0.0f, 2.0f);
            particles_.push_back(ccx+gv.x, ccy+gv.y, ccz+gv.z,
                                  cvx+rv.x, rv.y, cvz+rv.z,
                                  uniformRandom(0.8f, 1.2f));
        }
    }
}

void NBodySimulation::generatePlummerSphere(int n_particles, float scale_radius, float total_mass)
{
    particles_.clear();
    particles_.reserve(n_particles);

    if (total_mass < 0.0f) total_mass = static_cast<float>(n_particles);
    float particle_mass = total_mass / n_particles;

    float effective_radius = 2.4f * scale_radius;
    float sphere_volume = (4.f/3.f)*M_PI*effective_radius*effective_radius*effective_radius;
    float avg_spacing = std::cbrt(sphere_volume / n_particles);
    float softening = avg_spacing * 0.15f;
    params_.min_r2 = softening * softening;

    std::cout << "Generating Plummer Sphere:\n"
              << "  Particles: " << n_particles << "\n"
              << "  Scale radius: " << scale_radius << "\n"
              << "  Total mass: " << total_mass << "\n"
              << "  Particle mass: " << particle_mass << "\n"
              << "  Avg spacing: " << avg_spacing << "\n"
              << "  Softening: " << softening << " (min_r2 = " << params_.min_r2 << ")\n";
    
    for (int i = 0; i < n_particles; ++i) {
        float u = uniformRandom(0.001f, 1.0f);
        float r = scale_radius / std::sqrt(std::pow(u, -2.0f/3.0f) - 1.0f);
        auto pos = uniformRandomOnSphere(r);

        float x = 0.f, y = 0.f;
        while (true) {
            x = uniformRandom(0.f, 1.f);
            y = uniformRandom(0.f, 0.1f);
            float gx = x*x * std::pow(1.f - x*x, 3.5f);
            if (y < gx) break;
        }

        float r2a2 = r*r + scale_radius*scale_radius;
        float v_esc = std::sqrt(2.f * params_.G * total_mass / std::sqrt(r2a2));
        float v = x * v_esc;
        auto vel = uniformRandomOnSphere(v);

        particles_.push_back(pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, particle_mass);
    }

    // Remove bulk velocity
    float tot_px=0,tot_py=0,tot_pz=0,tot_m=0;
    for (size_t i = 0; i < particles_.size(); ++i) {
        tot_px += particles_.vel_x[i]*particles_.mass[i];
        tot_py += particles_.vel_y[i]*particles_.mass[i];
        tot_pz += particles_.vel_z[i]*particles_.mass[i];
        tot_m  += particles_.mass[i];
    }
    float bvx=tot_px/tot_m, bvy=tot_py/tot_m, bvz=tot_pz/tot_m;
    for (size_t i = 0; i < particles_.size(); ++i) {
        particles_.vel_x[i] -= bvx;
        particles_.vel_y[i] -= bvy;
        particles_.vel_z[i] -= bvz;
    }
    std::cout << "Plummer sphere generated successfully!\n";
}

void NBodySimulation::step()
{
    algorithm_->computeStep(particles_, params_);
    {
        std::lock_guard<std::mutex> lock(buffer_mutex_);
        render_buffer_ = particles_;
        buffer_ready_ = true;
    }
}

const ParticleSoA& NBodySimulation::getRenderBuffer()
{
    std::lock_guard<std::mutex> lock(buffer_mutex_);
    return render_buffer_;
}