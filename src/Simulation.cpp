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

void NBodySimulation::generateGalaxyDisk(int n_particles, float radius, float thickness)
{
    particles_.clear();
    particles_.reserve(n_particles);
    
    const float central_mass = 1000000.0f;
    particle_t center;
    center.position = {0, 0, 0};
    center.velocity = {0, 0, 0};
    center.mass = central_mass;
    particles_.push_back(center);
    
    for (int i = 1; i < n_particles; ++i) {
        particle_t p;
        
        float r = radius * std::sqrt(uniformRandom(0.0f, 1.0f));
        float theta = uniformRandom(0.0f, 2.0f * M_PI);
        float z = gaussianRandom(0.0f, thickness * 0.3f);
        
        p.position = {
            r * std::cos(theta),
            z,
            r * std::sin(theta)
        };
        
        float v_orbital = std::sqrt(params_.G * central_mass / r);
        float v_dispersion = 0.1f * v_orbital;
        
        p.velocity = {
            -v_orbital * std::sin(theta) + gaussianRandom(0.0f, v_dispersion),
            gaussianRandom(0.0f, v_dispersion * 0.5f),
            v_orbital * std::cos(theta) + gaussianRandom(0.0f, v_dispersion)
        };
        
        p.mass = uniformRandom(0.5f, 2.0f);
        particles_.push_back(p);
    }

    float disk_volume = M_PI * radius * radius * thickness;
    float volume_per_particle = disk_volume / n_particles;
    float avg_spacing = std::cbrt(volume_per_particle);
    
    float softening = avg_spacing * 0.1f;  // 10% of spacing
    params_.min_r2 = softening * softening;
    
    std::cout << "Disk: avg_spacing = " << avg_spacing 
              << ", softening = " << softening 
              << ", min_r2 = " << params_.min_r2 << "\n";
}

void NBodySimulation::generateClusters(int n_clusters, 
                                       int particles_per_cluster,
                                       float cluster_separation) 
{
    particles_.clear();
    particles_.reserve(n_clusters * particles_per_cluster);
    
    for (int c = 0; c < n_clusters; ++c) {
        float angle = (2.0f * M_PI * c) / n_clusters;
        vector3_t cluster_center = {
            cluster_separation * std::cos(angle),
            uniformRandom(-50.0f, 50.0f),
            cluster_separation * std::sin(angle)
        };
        
        float speed = 5.0f;
        vector3_t cluster_velocity = {
            -speed * std::cos(angle),
            0.0f,
            -speed * std::sin(angle)
        };
        
        for (int i = 0; i < particles_per_cluster; ++i) {
            particle_t p;
            
            float cluster_radius = 100.0f;
            p.position = cluster_center + gaussianRandomVec3(0.0f, cluster_radius);
            p.velocity = cluster_velocity + gaussianRandomVec3(0.0f, 2.0f);
            p.mass = uniformRandom(0.8f, 1.2f);
            
            particles_.push_back(p);
        }
    }
}

void NBodySimulation::generatePlummerSphere(int n_particles, 
                                           float scale_radius,
                                           float total_mass)
{
    particles_.clear();
    particles_.reserve(n_particles);
    
    // If total_mass not specified, use 1.0 per particle
    if (total_mass < 0.0f) {
        total_mass = static_cast<float>(n_particles);
    }
    float particle_mass = total_mass / n_particles;
    
    // Calculate softening based on Plummer sphere density
    // Effective radius containing ~90% of mass is about 2.4 * scale_radius
    float effective_radius = 2.4f * scale_radius;
    float sphere_volume = (4.0f/3.0f) * M_PI * 
                          effective_radius * effective_radius * effective_radius;
    float volume_per_particle = sphere_volume / n_particles;
    float avg_spacing = std::cbrt(volume_per_particle);
    float softening = avg_spacing * 0.15f;  // 15% for clusters
    
    params_.min_r2 = softening * softening;
    
    std::cout << "Generating Plummer Sphere:\n"
              << "  Particles: " << n_particles << "\n"
              << "  Scale radius: " << scale_radius << "\n"
              << "  Total mass: " << total_mass << "\n"
              << "  Particle mass: " << particle_mass << "\n"
              << "  Avg spacing: " << avg_spacing << "\n"
              << "  Softening: " << softening << " (min_r2 = " << params_.min_r2 << ")\n";
    
    // Generate particles
    for (int i = 0; i < n_particles; ++i)
    {
        particle_t p;
        p.mass = particle_mass;
        
        // ===== POSITION: Sample from Plummer density profile =====
        // Plummer density: ρ(r) = (3M / 4πa³) × (1 + r²/a²)^(-5/2)
        // Cumulative mass: M(r) = M × r³ / (r² + a²)^(3/2)
        // Inverse: r = a / sqrt(u^(-2/3) - 1), where u ~ Uniform(0,1)
        
        float u = uniformRandom(0.0f, 1.0f);
        float r = scale_radius / std::sqrt(std::pow(u, -2.0f/3.0f) - 1.0f);
        
        // Random direction on sphere
        p.position = uniformRandomOnSphere(r);
        
        // ===== VELOCITY: Sample from isotropic distribution function =====
        // For a self-consistent Plummer sphere in virial equilibrium:
        // Escape velocity: v_esc² = 2 × Φ(r), where Φ(r) = -GM / sqrt(r² + a²)
        // Distribution function requires rejection sampling
        
        float x = 0.0f;  // Ratio v / v_esc
        float y = 0.0f;
        
        // Rejection sampling for velocity magnitude
        // f(v) from Plummer's DF (Aarseth et al. 1974)
        while (true) {
            x = uniformRandom(0.0f, 1.0f);
            y = uniformRandom(0.0f, 0.1f);  // Upper bound on g(x)
            
            float g_x = x * x * std::pow(1.0f - x * x, 3.5f);
            
            if (y < g_x) {
                break;  // Accept this sample
            }
        }
        
        // Calculate escape velocity at this radius
        float r2_plus_a2 = r * r + scale_radius * scale_radius;
        float v_escape = std::sqrt(2.0f * params_.G * total_mass / std::sqrt(r2_plus_a2));
        
        // Actual velocity magnitude
        float v = x * v_escape;
        
        // Random direction for velocity (isotropic)
        p.velocity = uniformRandomOnSphere(v);
        
        particles_.push_back(p);
    }
    
    std::cout << "Plummer sphere generated successfully!\n";
}

void NBodySimulation::step()
{
    algorithm_->computeStep(particles_, params_);
    //std::cout << "computed step.\n";
    {
        std::lock_guard<std::mutex> lock(buffer_mutex_);
        render_buffer_ = particles_;  // or use std::swap
        buffer_ready_ = true;
    }
}

const std::vector<particle_t>& NBodySimulation::getRenderBuffer() 
{
    std::lock_guard<std::mutex> lock(buffer_mutex_);
    return render_buffer_;
}