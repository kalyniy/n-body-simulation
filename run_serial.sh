#!/bin/bash

# Configuration
RUNS_PER_CONFIG=5           # Number of runs per configuration
N_START=60000               # Starting particle count
N_END=100000                 # Ending particle count
STEPS=100                   # Number of simulation steps
EXECUTABLE="./bin/nbody_headless_seq"

# Generate N values: 10000, 20000, ..., 100000
n_values=()
for ((n=N_START; n<=N_END; n+=10000)); do
    n_values+=($n)
done

# Run the benchmarks
for n in "${n_values[@]}"; do
    echo "=== N=$n ==="
    for ((run=1; run<=RUNS_PER_CONFIG; run++)); do
        echo "    Run $run/$RUNS_PER_CONFIG"
        $EXECUTABLE --random $n --steps $STEPS --disable-output
    done
done