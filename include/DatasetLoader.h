#pragma once
#include <string>
#include <cmath>
#include <vector>
#include <filesystem>
#include <fstream>
#include "ParticleSoA.h"

class DatasetLoader
{
public:
    static size_t count_floats(const std::string &path);
    static void read_f32(const std::string &path, float *dst, size_t N);
    static void load_hacc_snapshot(ParticleSoA *particles, const std::string &dir);

    DatasetLoader() = delete;
};