#!/bin/bash

# Configuration
RUNS_PER_CONFIG=5           # Number of runs per configuration
N_START=100000               # Starting particle count
N_END=1000000                # Ending particle count
STEPS=100                   # Number of simulation steps
#EXECUTABLE="./bin/nbody_headless_seq"
EXECUTABLE="./bin/nbody_headless_cuda"
INCREMENT=100000
ALGORITHM="bh"

# Generate N values: 10000, 20000, ..., 100000
n_values=()
for ((n=N_START; n<=N_END; n+=INCREMENT)); do
    n_values+=($n)
done

# Run the benchmarks
for n in "${n_values[@]}"; do
    echo "=== N=$n ==="
    for ((run=1; run<=RUNS_PER_CONFIG; run++)); do
        echo "    Run $run/$RUNS_PER_CONFIG"
        $EXECUTABLE --algo $ALGORITHM --plummer $n --steps $STEPS --disable-output
    done
done