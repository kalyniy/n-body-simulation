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
        fp_ = std::fopen(path.c_str(), "w");
        if (fp_)
        {
            switch (mode)
            {
            case RunMode::Run:
                printRunHeader();
                break;
            case RunMode::Step:
                printStepHeader();
                break;
            default:
                break;
            }
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

        return (29LL * (n * n) - 5LL * n) / 2;
    }

    void log(size_t n_particles, size_t n_steps, double total_time) /*size_t N, double update_sec, double gflops = 0.0, const char *note = ""*/
    {
        constexpr double giga = 1e9; // Use decimal giga (10^9)

        long long operations_count = calculate_total_number_of_flops_per_step(n_particles) * n_steps;
        printf("Total operations: %lld\n", operations_count);
        printf("Total time: %lf seconds\n", total_time);

        double flops = (double)operations_count / total_time;
        printf("FLOP/s: %.2e\n", flops);

        double gflops = flops / giga;
        printf("GFLOP/s: %.2lf\n", gflops);

        if (fp_)
        {
            std::fprintf(fp_, "%12zu %12zu %12.6lf %14.2lf\n",
                         n_particles, n_steps, total_time, gflops);
        }
    }

    bool ok() const { return fp_ != nullptr; }

private:
    std::FILE *fp_ = nullptr;

    void printStepHeader()
    {
        std::fprintf(fp_, "N-Body simulation Performance Per Step\n");
        std::fprintf(fp_, "============================\n\n");
        std::fprintf(fp_, "%6s %12s %14s %15s\n", "N Particles", "", "GFLOPS", "Notes");
        std::fprintf(fp_, "---------------------------------------------------------------------\n");
    }

    void printRunHeader()
    {
        std::fprintf(fp_, "N-Body simulation Performance Per Run\n");
        std::fprintf(fp_, "============================\n\n");
        std::fprintf(fp_, "%12s %12s %12s %14s\n", "N Particles", "N Steps", "Time (s)", "GFLOPS");
        std::fprintf(fp_, "---------------------------------------------------------------------\n");
    }
};
