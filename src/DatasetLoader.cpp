#include "DatasetLoader.h"
#include <cstdio>
#include <memory>

size_t DatasetLoader::count_floats(const std::string &path) {
    namespace fs = std::filesystem;
    auto bytes = fs::file_size(path);
    if (bytes % sizeof(float)) throw std::runtime_error("Not multiple of 4B: " + path);
    return bytes / sizeof(float);
}

void DatasetLoader::read_f32(const std::string &path, float *dst, size_t N) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + path);
    f.read(reinterpret_cast<char*>(dst), N * sizeof(float));
    if (f.gcount() != static_cast<std::streamsize>(N * sizeof(float)))
        throw std::runtime_error("Short read on " + path);
}

void DatasetLoader::load_hacc_snapshot(ParticleSoA *particles, const std::string &dir) {
    const std::string fx=dir+"/xx.f32",fy=dir+"/yy.f32",fz=dir+"/zz.f32";
    const std::string fvx=dir+"/vx.f32",fvy=dir+"/vy.f32",fvz=dir+"/vz.f32";
    size_t N = count_floats(fx);
    std::printf("[HACC] Found %zu particles.\n", N);

    particles->resize(N);
    read_f32(fx, particles->pos_x.data(), N);
    read_f32(fy, particles->pos_y.data(), N);
    read_f32(fz, particles->pos_z.data(), N);
    read_f32(fvx, particles->vel_x.data(), N);
    read_f32(fvy, particles->vel_y.data(), N);
    read_f32(fvz, particles->vel_z.data(), N);
    // Mass defaults to 1.0 for HACC
    std::fill(particles->mass.begin(), particles->mass.end(), 1.0f);
    particles->zeroAccelerations();
}