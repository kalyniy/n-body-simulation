#pragma once
#include <vector>
#include <cstring>
#include <cassert>

struct ParticleSoA
{
    std::vector<float> pos_x, pos_y, pos_z;
    std::vector<float> vel_x, vel_y, vel_z;
    std::vector<float> acc_x, acc_y, acc_z;
    std::vector<float> mass;

    size_t size() const { return mass.size(); }
    bool   empty() const { return mass.empty(); }

    void resize(size_t n)
    {
        pos_x.resize(n); pos_y.resize(n); pos_z.resize(n);
        vel_x.resize(n); vel_y.resize(n); vel_z.resize(n);
        acc_x.resize(n); acc_y.resize(n); acc_z.resize(n);
        mass.resize(n);
    }

    void reserve(size_t n)
    {
        pos_x.reserve(n); pos_y.reserve(n); pos_z.reserve(n);
        vel_x.reserve(n); vel_y.reserve(n); vel_z.reserve(n);
        acc_x.reserve(n); acc_y.reserve(n); acc_z.reserve(n);
        mass.reserve(n);
    }

    void clear()
    {
        pos_x.clear(); pos_y.clear(); pos_z.clear();
        vel_x.clear(); vel_y.clear(); vel_z.clear();
        acc_x.clear(); acc_y.clear(); acc_z.clear();
        mass.clear();
    }

    void push_back(float px, float py, float pz,
                    float vx, float vy, float vz,
                    float m)
    {
        pos_x.push_back(px); pos_y.push_back(py); pos_z.push_back(pz);
        vel_x.push_back(vx); vel_y.push_back(vy); vel_z.push_back(vz);
        acc_x.push_back(0.f); acc_y.push_back(0.f); acc_z.push_back(0.f);
        mass.push_back(m);
    }

    void zeroAccelerations()
    {
        std::memset(acc_x.data(), 0, acc_x.size() * sizeof(float));
        std::memset(acc_y.data(), 0, acc_y.size() * sizeof(float));
        std::memset(acc_z.data(), 0, acc_z.size() * sizeof(float));
    }
};