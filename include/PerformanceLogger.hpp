#pragma once
#include <cstdio>
#include <string>

class PerformanceLogger
{
public:
    explicit PerformanceLogger(const std::string &path)
    {
        fp_ = std::fopen(path.c_str(), "w");
        if (fp_)
        {
            std::fprintf(fp_, "N-Body simulation Performance\n");
            std::fprintf(fp_, "============================\n\n");
            std::fprintf(fp_, "%6s %12s %14s %15s\n", "N", "Update (s)", "GFLOPS", "Notes");
            std::fprintf(fp_, "---------------------------------------------------------------------\n");
        }
    }
    ~PerformanceLogger()
    {
        if (fp_)
            std::fclose(fp_);
    }
    
    long long calculate_total_number_of_flops_per_step(size_t n_particles)
    {
        long long n = (long long)n_particles;
        
        return (29LL*(n*n) - 5LL*n)/2;
    }

    void log(size_t n_particles, size_t n_steps, double total_time) /*size_t N, double update_sec, double gflops = 0.0, const char *note = ""*/
    {
        constexpr size_t giga = 1024*1024*1024;

        long long operations_count = calculate_total_number_of_flops_per_step(n_particles);
        printf("Total operations per step: %lld\n", operations_count);;

        double flops = (double)operations_count/total_time;

        auto gflops = flops / giga;

        if (fp_)
        {
            std::fprintf(fp_, "%6zu %6zu %12.6lf %14.2lf\n", n_particles, n_steps, total_time, gflops);
        }
    }

    bool ok() const { return fp_ != nullptr; }

private:
    std::FILE *fp_ = nullptr;
};
