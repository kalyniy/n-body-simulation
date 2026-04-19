// CudaMisc.h
#pragma once
#include <stdio.h>
#include <cuda_runtime.h>

static inline void cudaCheck(cudaError_t e, const char* msg) {
    if (e != cudaSuccess) {
        fprintf(stderr, "CUDA error: %s: %s\n", msg, cudaGetErrorString(e));
        #ifdef USE_MPI
        MPI_Finalize();
        #endif
        std::exit(EXIT_FAILURE);
    }
}