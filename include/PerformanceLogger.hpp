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

    void log(size_t N, double update_sec, double gflops = 0.0, const char *note = "")
    {
        if (fp_)
            std::fprintf(fp_, "%6zu %12.6f %14.2f %15s\n", N, update_sec, gflops, note);
    }

    bool ok() const { return fp_ != nullptr; }

private:
    std::FILE *fp_ = nullptr;
};
