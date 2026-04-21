# N-Body Simulation
 
A high-performance gravitational N-body simulation with sequential, MPI, CUDA, and hybrid CUDA+MPI backends, using either a naïve O(N²) or Barnes-Hut O(N log N) algorithm. Simulation output is written to a binary checkpoint file and visualized in ParaView via a VTK converter.
 
---
 
## Table of Contents
 
1. [Requirements](#requirements)
2. [Project Structure](#project-structure)
3. [Building](#building)
4. [Running](#running)
5. [Visualization with ParaView](#visualization-with-paraview)
6. [Initial Conditions](#initial-conditions)
7. [Physics Parameters](#physics-parameters)
8. [Performance Logging](#performance-logging)
9. [Quick Reference](#quick-reference)
---
 
## Requirements
 
| Dependency | Required for |
|---|---|
| `g++` with C++20 support | All builds |
| `mpic++` / OpenMPI | `mpi`, `cuda_mpi` builds |
| CUDA Toolkit + `nvcc` | `cuda`, `cuda_mpi` builds |
| ParaView | Visualization |
 
---
 
## Project Structure
 
```
n-body-simulation/
├── Makefile                    # Main build file
├── include/                    # Header files
├── src/
│   ├── apps/
│   │   └── headless_main.cpp   # Main entry point
│   └── tools/
│       └── checkpoint_to_vtk.cpp
├── bin/                        # Compiled binaries (generated)
│   ├── nbody_headless_seq
│   ├── nbody_headless_mpi
│   ├── nbody_headless_cuda
│   ├── nbody_headless_cuda_mpi
│   └── checkpoint_to_vtk
├── data/                       # Performance CSV results
└── simulation_output.bin       # Binary checkpoint (generated)
```
 
---
 
## Building
 
The Makefile supports four build modes selected via `MODE=`.
 
```bash
make            # Sequential CPU (default, same as MODE=seq)
make MODE=mpi        # MPI distributed CPU
make MODE=cuda       # Single-GPU CUDA
make MODE=cuda_mpi   # Multi-GPU CUDA + MPI hybrid
make clean           # Remove all build artifacts
```
 
To override the CUDA architecture (default: auto-detected via `nvidia-smi`, fallback `sm_70`):
 
```bash
make MODE=cuda CUDA_ARCH=86   # e.g., for Ampere RTX 30xx
```
 
To build only the VTK converter tool:
 
```bash
make vtk
```
 
---
 
## Running
 
### Sequential
 
```bash
./bin/nbody_headless_seq --algo naive --random 10000 --steps 500 --dt 0.1 --out output.txt
```
 
### MPI
 
```bash
mpirun -np 4 ./bin/nbody_headless_mpi --algo naive --random 10000 --steps 500
```
 
Or use the provided script:
 
```bash
bash run_mpi.sh
```
 
### CUDA
 
```bash
./bin/nbody_headless_cuda --algo bh --plummer 50000 --steps 1000
```
 
### CUDA + MPI
 
```bash
bash run_cuda_mpi.sh
```
 
---
 
### Command-Line Options
 
| Flag | Description | Default |
|---|---|---|
| `--algo <naive\|bh>` | Algorithm: naïve O(N²) or Barnes-Hut O(N log N) | `naive` |
| `--theta <float>` | Barnes-Hut opening angle (smaller = more accurate) | `0.5` |
| `--steps <N>` | Number of simulation steps | `1000` |
| `--dt <float>` | Time step size | `0.5` |
| `--G <float>` | Gravitational constant | `1.0` |
| `--softening <float>` | Gravitational softening length | `2.0` |
| `--out <file>` | Performance log output file | `output.txt` |
| `--log-every <N>` | Print progress every N steps | `100` |
| `--disable-output` | Skip writing checkpoint binary | — |
 
Initial condition flags are described in the [Initial Conditions](#initial-conditions) section below.
 
---
 
## Visualization with ParaView
 
Simulation state is written to `simulation_output.bin` after each run. Convert it to a VTK time series with:
 
```bash
./bin/checkpoint_to_vtk simulation_output.bin ./vtk_output --every 1
```
 
| Option | Description |
|---|---|
| `--every N` | Export every Nth step (useful for large runs) |
| `--start N` | Start from step N |
| `--end N` | Stop at step N |
 
This produces one `.vtk` file per exported step plus a `nbody.vtk.series` index file inside `./vtk_output/`.
 
### Loading in ParaView
 
1. Open ParaView.
2. **File → Open** → select `vtk_output/nbody.vtk.series`.
3. Click **Apply** in the Properties panel.
4. In the toolbar, change **Coloring** from `Solid Color` to `mass`.
5. Press **Play** to animate the time series.
> **Tip:** Increase *Point Size* in the Properties panel if particles appear too small.
 
---
 
## Initial Conditions
 
Exactly one initial condition flag must be provided per run.
 
| Flag | Description |
|---|---|
| `--random <N>` | N particles uniformly distributed in the domain |
| `--plummer <N>` | Plummer sphere (self-gravitating equilibrium cluster) |
| `--galaxy <N>` | Rotating disk galaxy with central mass |
| `--clusters <N>` | 3 colliding clusters, N particles each |
| `--solar` | 9-body solar system analog |
| `--hacc <dir>` | Load a HACC cosmological snapshot from directory |
 
### Examples
 
```bash
# 50,000-particle Plummer sphere, Barnes-Hut, 2000 steps
./bin/nbody_headless_seq --algo bh --plummer 50000 --steps 2000
 
# Galaxy disk, CUDA, log every 50 steps
./bin/nbody_headless_cuda --algo bh --galaxy 100000 --steps 1000 --log-every 50
 
# Three colliding clusters, MPI across 8 ranks
mpirun -np 8 ./bin/nbody_headless_mpi --algo naive --clusters 5000 --steps 500
 
# Solar system, sequential, no checkpoint output
./bin/nbody_headless_seq --solar --steps 10000 --disable-output
```
 
---
 
## Physics Parameters
 
The softening length (`--softening`) is automatically set to a physically sensible value when using `--plummer` or `--galaxy`. For `--random` and other modes you may want to tune it manually.
 
Barnes-Hut theta (`--theta`) controls the accuracy–speed trade-off:
 
- `0.0` → exact (equivalent to naïve)
- `0.5` → recommended default
- `1.0` → fast but less accurate
---
 
## Performance Logging
 
Each run appends a row to the file specified by `--out` (default `output.txt`), recording particle count, step count, wall-clock time, and process count. Pre-collected results are in `data/`:
 
| File | Description |
|---|---|
| `seq-naive.csv` | Sequential naïve |
| `seq-barnes-hut.csv` | Sequential Barnes-Hut |
| `mpi-naive.csv` | MPI naïve |
| `mpi-barnes-hut.csv` | MPI Barnes-Hut |
| `seq-cuda-naive.csv` | CUDA naïve |
| `seq-cuda-barnes-hut.csv` | CUDA Barnes-Hut |
| `seq-cuda-barnes-hut-v2.csv` | CUDA Barnes-Hut (revised) |
 
---
 
## Quick Reference
 
```bash
# Build everything
make seq && make mpi && make cuda
 
# Run a quick test (sequential, 1000 particles, 100 steps)
./bin/nbody_headless_seq --algo naive --random 1000 --steps 100
 
# Convert last run to VTK and open in ParaView
./bin/checkpoint_to_vtk simulation_output.bin ./vtk_output --every 1
 
# Get help
./bin/nbody_headless_seq --help
```