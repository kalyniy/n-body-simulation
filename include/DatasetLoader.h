#pragma once
#include <string>
#include <bits/stdc++.h>
#include <vector>
#include <filesystem>

#include "Particle.hpp"

class DatasetLoader
{
public:
    static size_t count_floats(const std::string &path);
    static void read_f32(const std::string &path, float *dst, size_t N);
    static void load_hacc_snapshot(std::vector<particle_t> *particles, const std::string &dir);

    DatasetLoader() = delete;
};