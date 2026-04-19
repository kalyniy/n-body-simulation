#pragma once
#include <iostream>
#include <string>
#include <mutex>

#include "SimulationOutput.h"
#include "ParticleSoA.h"

class CheckpointManager {
private:
    std::string filePath;
    static CheckpointManager* instancePtr;
    static std::mutex mtx;
    CheckpointManager() {}
public:
    CheckpointManager(const CheckpointManager&) = delete;
    static CheckpointManager* getInstance() {
        if (!instancePtr) {
            std::lock_guard<std::mutex> lock(mtx);
            if (!instancePtr) instancePtr = new CheckpointManager();
        }
        return instancePtr;
    }
    void setFilePath(const std::string& path) { 
        filePath = path;
    }
    void write_header(SimulationOutputHeader header);
    void write_masses(const ParticleSoA& particles);
    void increment_passed_steps();
    void write_step(const ParticleSoA& particles);

    SimulationOutputHeader read_header();
    void read_masses(float* masses_out, size_t n_particles);
    size_t read_step(float* positions_out, size_t step_index);
};