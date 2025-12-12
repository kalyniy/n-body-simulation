#!/bin/bash

# Configuration
MAX_NP=8                    # Maximum number of processes (change to 16, etc.)
RUNS_PER_CONFIG=5           # Number of runs per configuration
N_START=70000               # Starting particle count
N_END=100000                # Ending particle count
STEPS=100                   # Number of simulation steps
EXECUTABLE="./bin/nbody_headless_mpi"

# Generate np values: 1, 2, 4, 6, 8, ...
np_values=(1)
for ((np=2; np<=MAX_NP; np+=2)); do
    np_values+=($np)
done

# Generate N values: 10000, 20000, ..., 100000
n_values=()
for ((n=N_START; n<=N_END; n+=10000)); do
    n_values+=($n)
done

# Run the benchmarks
for n in "${n_values[@]}"; do
    echo "=== N=$n ==="
    for np in "${np_values[@]}"; do
        echo "  np=$np"
        for ((run=1; run<=RUNS_PER_CONFIG; run++)); do
            echo "    Run $run/$RUNS_PER_CONFIG"
            mpirun -np $np $EXECUTABLE --random $n --steps $STEPS --disable-output
        done
    done
done