#include "DatasetLoader.h"
#include <filesystem>

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

ParticlesSOA DatasetLoader::load_hacc_snapshot(const std::string &dir)
{
    const std::string fx = dir + "/xx.f32", fy = dir + "/yy.f32", fz = dir + "/zz.f32";
    const std::string fvx = dir + "/vx.f32", fvy = dir + "/vy.f32", fvz = dir + "/vz.f32";

    size_t Nx = count_floats(fx), Ny = count_floats(fy), Nz = count_floats(fz);
    size_t Nvx = count_floats(fvx), Nvy = count_floats(fvy), Nvz = count_floats(fvz);
    if (!(Nx == Ny && Nx == Nz && Nx == Nvx && Nx == Nvy && Nx == Nvz))
        throw std::runtime_error("Array sizes differ");

    ParticlesSOA P;
    P.N = Nx;
    P.x = std::make_unique<float[]>(P.N);
    P.y = std::make_unique<float[]>(P.N);
    P.z = std::make_unique<float[]>(P.N);
    P.vx = std::make_unique<float[]>(P.N);
    P.vy = std::make_unique<float[]>(P.N);
    P.vz = std::make_unique<float[]>(P.N);

    read_f32(fx, P.x.get(), P.N);
    read_f32(fy, P.y.get(), P.N);
    read_f32(fz, P.z.get(), P.N);
    read_f32(fvx, P.vx.get(), P.N);
    read_f32(fvy, P.vy.get(), P.N);
    read_f32(fvz, P.vz.get(), P.N);
    return P;
}