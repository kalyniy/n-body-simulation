#pragma once
#include <cstdio>
#include <string>

enum RunMode
{
    Run,
    Step
};

class PerformanceLogger
{
public:
    explicit PerformanceLogger(const std::string &path, RunMode mode)
    {
        fp_ = std::fopen(path.c_str(), "a");
    }

    ~PerformanceLogger()
    {
        if (fp_)
            std::fclose(fp_);
    }

    void log(size_t n_particles, size_t n_steps, double total_time, size_t workers_count = 1)
    {
        constexpr double giga = 1e9; // Use decimal giga (10^9)

        if (fp_)
        {
            std::fprintf(fp_, "%zu,%zu,%lf,%zu\n", n_particles, n_steps, total_time, workers_count);
        }
    }

    bool ok() const { return fp_ != nullptr; }

private:
    std::FILE *fp_ = nullptr;
};
