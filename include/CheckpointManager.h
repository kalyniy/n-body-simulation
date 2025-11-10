#pragma once
#include <iostream>
#include <string>
#include <mutex>

#include "SimulationOutput.h"
#include "Particle.hpp"

class CheckpointManager {
private:
    std::string filePath;

    // Static pointer to the Singleton instance
    static CheckpointManager* instancePtr;

    // Mutex to ensure thread safety
    static std::mutex mtx;

    // Private Constructor
    CheckpointManager() {}

public:
    // Deleting the copy constructor to prevent copies
    CheckpointManager(const CheckpointManager& obj) = delete;

    // Static method to get the Singleton instance
    static CheckpointManager* getInstance() {
        if (instancePtr == nullptr) {
            std::lock_guard<std::mutex> lock(mtx);
            if (instancePtr == nullptr) {
                instancePtr = new CheckpointManager();
            }
        }
        return instancePtr;
    }

    // Method to set values
    void setFilePath(const std::string& name) {
        this->filePath = filePath;
    }

    void write_header(SimulationOutputHeader header);
    void increment_passed_steps();
    void write_step(particle_t* particles, int count);

    SimulationOutputHeader read_header();
    void read_step();
};
