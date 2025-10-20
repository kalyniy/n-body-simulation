#pragma once
#include <random>
#include <cmath>
#include "Particle.hpp"

// Thread-local random engine for better performance in parallel code
inline std::mt19937& getRandomEngine() {
    static thread_local std::mt19937 engine(std::random_device{}());
    return engine;
}

// Uniform random float in range [min, max)
inline float uniformRandom(float min, float max) {
    std::uniform_real_distribution<float> dist(min, max);
    return dist(getRandomEngine());
}

// Gaussian (normal) random with mean and standard deviation
inline float gaussianRandom(float mean, float stddev) {
    std::normal_distribution<float> dist(mean, stddev);
    return dist(getRandomEngine());
}

// Gaussian random vector (for 3D positions/velocities)
inline vector3_t gaussianRandomVec3(float mean, float stddev) {
    return {
        gaussianRandom(mean, stddev),
        gaussianRandom(mean, stddev),
        gaussianRandom(mean, stddev)
    };
}

// Uniform random vector in a box
inline vector3_t uniformRandomVec3(float min_x, float max_x,
                                   float min_y, float max_y,
                                   float min_z, float max_z) {
    return {
        uniformRandom(min_x, max_x),
        uniformRandom(min_y, max_y),
        uniformRandom(min_z, max_z)
    };
}

// Uniform random vector in a cube (convenience)
inline vector3_t uniformRandomVec3(float min, float max) {
    return uniformRandomVec3(min, max, min, max, min, max);
}

// Uniform random point on a sphere surface (radius r)
inline vector3_t uniformRandomOnSphere(float radius) {
    // Marsaglia method
    float x, y, s;
    do {
        x = uniformRandom(-1.0f, 1.0f);
        y = uniformRandom(-1.0f, 1.0f);
        s = x*x + y*y;
    } while (s >= 1.0f);
    
    float t = 2.0f * std::sqrt(1.0f - s);
    return {
        radius * x * t,
        radius * y * t,
        radius * (1.0f - 2.0f * s)
    };
}

// Uniform random point inside a sphere (radius r)
inline vector3_t uniformRandomInSphere(float radius) {
    // Use cube root for uniform volume distribution
    float r = radius * std::cbrt(uniformRandom(0.0f, 1.0f));
    return uniformRandomOnSphere(r);
}

// Uniform random point in a disk (radius r, on XZ plane)
inline vector3_t uniformRandomInDisk(float radius, float thickness = 0.0f) {
    float r = radius * std::sqrt(uniformRandom(0.0f, 1.0f));
    float theta = uniformRandom(0.0f, 2.0f * M_PI);
    float z = (thickness > 0.0f) ? uniformRandom(-thickness, thickness) : 0.0f;
    
    return {
        r * std::cos(theta),
        z,
        r * std::sin(theta)
    };
}