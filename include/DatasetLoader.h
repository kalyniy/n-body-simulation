#pragma once
#include <string>
#include <bits/stdc++.h>

struct ParticlesSOA
{
    size_t N{};
    std::unique_ptr<float[]> x, y, z, vx, vy, vz;
};

class DatasetLoader
{
public:
    static size_t count_floats(const std::string &path);
    static void read_f32(const std::string &path, float *dst, size_t N);
    static ParticlesSOA load_hacc_snapshot(const std::string &dir);

    DatasetLoader() = delete;
};