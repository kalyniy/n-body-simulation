#include "DatasetLoader.h"

#include <stdio.h>

using std::size_t;

size_t DatasetLoader::count_floats(const std::string &path)
{
    namespace fs = std::filesystem;
    auto bytes = fs::file_size(path);
    if (bytes % sizeof(float))
        throw std::runtime_error("File not multiple of 4B: " + path);
    return bytes / sizeof(float);
}

void DatasetLoader::read_f32(const std::string &path, float *dst, size_t N)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        throw std::runtime_error("Cannot open " + path);
    f.read(reinterpret_cast<char *>(dst), N * sizeof(float));
    if (f.gcount() != static_cast<std::streamsize>(N * sizeof(float)))
        throw std::runtime_error("Short read on " + path);
}

void DatasetLoader::load_hacc_snapshot(std::vector<particle_t> *particles, const std::string &dir)
{
    const std::string fx = dir + "/xx.f32", fy = dir + "/yy.f32", fz = dir + "/zz.f32";
    const std::string fvx = dir + "/vx.f32", fvy = dir + "/vy.f32", fvz = dir + "/vz.f32";

    size_t Nx = count_floats(fx), Ny = count_floats(fy), Nz = count_floats(fz);
    size_t Nvx = count_floats(fvx), Nvy = count_floats(fvy), Nvz = count_floats(fvz);
    if (!(Nx == Ny && Nx == Nz && Nx == Nvx && Nx == Nvy && Nx == Nvz))
        throw std::runtime_error("Array sizes differ");

    printf("[HACC Dataset] Found %lu particles.\n", Nx);

    particles->reserve(Nx);

    auto vectorX = std::make_unique<float[]>(Nx);
    auto vectorY = std::make_unique<float[]>(Nx);
    auto vectorZ = std::make_unique<float[]>(Nx);
    auto vectorVx = std::make_unique<float[]>(Nx);
    auto vectorVy = std::make_unique<float[]>(Nx);
    auto vectorVz = std::make_unique<float[]>(Nx);

    printf("reading vector of x...\n");
    read_f32(fx, vectorX.get(), Nx);

    printf("reading vector of y...\n");
    read_f32(fy, vectorY.get(), Nx);

    printf("reading vector of z...\n");
    read_f32(fz, vectorZ.get(), Nx);

    printf("reading vector of vx...\n");
    read_f32(fvx, vectorVx.get(), Nx);

    printf("reading vector of vy...\n");
    read_f32(fvy, vectorVy.get(), Nx);

    printf("reading vector of vz...\n");
    read_f32(fvz, vectorVz.get(), Nx);

    for (size_t i = 0; i < Nx; i++)
    {
        particle_t p = {0};
        p.position.x = vectorX[i];
        p.position.y = vectorY[i];
        p.position.z = vectorZ[i];

        particles->push_back(p);
    }

    vectorX.release();
    vectorY.release();
    vectorZ.release();
    vectorVx.release();
    vectorVy.release();
    vectorVz.release();
}