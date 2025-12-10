// checkpoint_to_vtk.cpp
// Converts n-body simulation checkpoint files to VTK format for ParaView
// 
// Build:   g++ -std=c++17 -O2 -o checkpoint_to_vtk checkpoint_to_vtk.cpp
// Usage:   ./checkpoint_to_vtk simulation_output.bin ./vtk_output
//
// This will create vtk_output/ directory with:
//   - nbody_000000.vtk, nbody_000001.vtk, ... (one per timestep)
//   - nbody.vtk.series (for easy loading in ParaView)

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <iomanip>
#include <sstream>
#include <cstdint>

#include "SimulationOutput.h"

void writeVTK(const std::string& filename,
              const float* positions,
              const float* masses,
              size_t n_particles) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open " << filename << " for writing\n";
        return;
    }

    // VTK Legacy ASCII header
    file << "# vtk DataFile Version 3.0\n";
    file << "N-Body Simulation\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    // Points
    file << "POINTS " << n_particles << " float\n";
    file << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < n_particles; i++) {
        file << positions[i * 3 + 0] << " "
             << positions[i * 3 + 1] << " "
             << positions[i * 3 + 2] << "\n";
    }

    // Vertices (required for point visualization in ParaView)
    file << "VERTICES " << n_particles << " " << (n_particles * 2) << "\n";
    for (size_t i = 0; i < n_particles; i++) {
        file << "1 " << i << "\n";
    }

    // Point data - mass for coloring
    file << "POINT_DATA " << n_particles << "\n";
    file << "SCALARS mass float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < n_particles; i++) {
        file << masses[i] << "\n";
    }

    file.close();
}

std::string makeFilename(const std::string& outputDir, size_t step, int padding = 6) {
    std::ostringstream ss;
    ss << outputDir << "/nbody_" << std::setfill('0') << std::setw(padding) << step << ".vtk";
    return ss.str();
}

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " <checkpoint_file> <output_dir> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  --every N     Export every Nth step (default: 1)\n";
    std::cout << "  --start N     Start from step N (default: 0)\n";
    std::cout << "  --end N       End at step N (default: all)\n";
    std::cout << "  --help        Show this help message\n\n";
    std::cout << "Example:\n";
    std::cout << "  " << progName << " simulation_output.bin ./vtk_output --every 10\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }

    std::string checkpointFile = argv[1];
    std::string outputDir = argv[2];
    
    size_t exportEvery = 1;
    size_t startStep = 0;
    size_t endStep = SIZE_MAX;

    // Parse optional arguments
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--every" && i + 1 < argc) {
            exportEvery = std::stoul(argv[++i]);
        } else if (arg == "--start" && i + 1 < argc) {
            startStep = std::stoul(argv[++i]);
        } else if (arg == "--end" && i + 1 < argc) {
            endStep = std::stoul(argv[++i]);
        } else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
    }

    // Create output directory
    mkdir(outputDir.c_str(), 0755);

    // Open checkpoint file
    std::ifstream file(checkpointFile, std::ios::in | std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot open checkpoint file: " << checkpointFile << "\n";
        return 1;
    }

    // Read header
    SimulationOutputHeader header;
    file.read(reinterpret_cast<char*>(&header), sizeof(SimulationOutputHeader));
    
    std::cout << "Checkpoint Info:\n";
    std::cout << "  Particles:       " << header.n_particles << "\n";
    std::cout << "  Target steps:    " << header.target_steps << "\n";
    std::cout << "  Completed steps: " << header.passed_steps << "\n\n";

    // Read masses (stored once after header)
    std::vector<float> masses(header.n_particles);
    file.read(reinterpret_cast<char*>(masses.data()), header.n_particles * sizeof(float));
    
    std::cout << "Loaded masses (first 5): ";
    for (size_t i = 0; i < std::min<size_t>(5, header.n_particles); i++) {
        std::cout << masses[i] << " ";
    }
    std::cout << "\n\n";

    // Clamp endStep
    if (endStep > header.passed_steps) {
        endStep = header.passed_steps;
    }

    // Allocate position buffer
    std::vector<float> positions(header.n_particles * 3);
    size_t step_size = header.n_particles * 3 * sizeof(float);

    // Track exported steps for series file
    std::vector<std::pair<std::string, size_t>> exportedSteps;

    // Convert each step
    std::cout << "Converting steps to VTK (every " << exportEvery << " steps)...\n";
    
    size_t exportCount = 0;
    for (size_t step = 0; step < header.passed_steps; step++) {
        // Read positions for this step (we must read sequentially)
        file.read(reinterpret_cast<char*>(positions.data()), step_size);
        
        if (!file) {
            std::cerr << "Error reading step " << step << "\n";
            break;
        }

        // Skip if not in range or not on interval
        if (step < startStep || step >= endStep || step % exportEvery != 0) {
            continue;
        }

        // Write VTK file
        std::string vtkFilename = "nbody_" + 
            std::string(6 - std::to_string(step).length(), '0') + 
            std::to_string(step) + ".vtk";
        std::string vtkPath = outputDir + "/" + vtkFilename;
        
        writeVTK(vtkPath, positions.data(), masses.data(), header.n_particles);
        exportedSteps.push_back({vtkFilename, step});
        exportCount++;
        
        if (exportCount % 100 == 0 || step == endStep - 1) {
            std::cout << "  Exported step " << step << " (" << exportCount << " files)\n";
        }
    }

    file.close();

    // Create ParaView series file for easy loading
    std::string seriesFile = outputDir + "/nbody.vtk.series";
    std::ofstream series(seriesFile);
    series << "{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n";
    for (size_t i = 0; i < exportedSteps.size(); i++) {
        series << "    { \"name\" : \"" << exportedSteps[i].first 
               << "\", \"time\" : " << exportedSteps[i].second << " }";
        if (i < exportedSteps.size() - 1) series << ",";
        series << "\n";
    }
    series << "  ]\n}\n";
    series.close();

    std::cout << "\n==========================================\n";
    std::cout << "Done! Exported " << exportCount << " VTK files to: " << outputDir << "/\n";
    std::cout << "==========================================\n\n";
    
    std::cout << "Next steps in ParaView:\n";
    std::cout << "  1. Open ParaView\n";
    std::cout << "  2. File -> Open -> select: " << outputDir << "/nbody.vtk.series\n";
    std::cout << "  3. Click 'Apply' in the Properties panel (left side)\n";
    std::cout << "  4. Change 'Coloring' dropdown from 'Solid Color' to 'mass'\n";
    std::cout << "  5. Click play button to animate\n\n";
    std::cout << "Tip: Increase 'Point Size' in Properties panel if particles are too small\n";

    return 0;
}