# ===========================================================================
#  N-Body Simulation — Makefile
#
#  Supported modes (set via MODE= on the command line):
#      seq        — sequential CPU                       (default)
#      mpi        — MPI-distributed CPU
#      cuda       — single-GPU CUDA
#      cuda_mpi   — multi-GPU CUDA + MPI hybrid
#
#  Examples:
#      make                          # builds seq
#      make MODE=mpi                 # builds MPI
#      make MODE=cuda                # builds single-GPU CUDA
#      make MODE=cuda_mpi            # builds CUDA + MPI hybrid
#      make clean
# ===========================================================================

MODE ?= seq

MPI_CFLAGS  :=
MPI_LDFLAGS :=

# ---- Toolchain selection per mode ----------------------------------------

ifeq ($(MODE),cuda_mpi)
    CXX       = mpic++
    NVCC      = nvcc
    MPI_CFLAGS  := $(shell mpic++ --showme:compile)
    MPI_LDFLAGS := $(shell mpic++ --showme:link)
    CXXFLAGS  = -std=c++20 -Wall -O3 -I./include -DUSE_MPI -DUSE_CUDA
    NVCCFLAGS = -std=c++17 -O3 -I./include -DUSE_MPI -DUSE_CUDA -Xcompiler -Wall $(MPI_CFLAGS)
    MODE_SUFFIX = _cuda_mpi
    OBJ_DIR     = bin/cuda_mpi
    USE_CUDA_FLAG = 1
    USE_MPI_FLAG  = 1
else ifeq ($(MODE),cuda)
    CXX       = g++
    NVCC      = nvcc
    CXXFLAGS  = -std=c++20 -Wall -O3 -I./include -DUSE_CUDA
    NVCCFLAGS = -std=c++17 -O3 -I./include -DUSE_CUDA -Xcompiler -Wall
    MODE_SUFFIX = _cuda
    OBJ_DIR     = bin/cuda
    USE_CUDA_FLAG = 1
    USE_MPI_FLAG  = 0
else ifeq ($(MODE),mpi)
    CXX       = mpic++
    CXXFLAGS  = -std=c++20 -Wall -O3 -I./include -DUSE_MPI
    MODE_SUFFIX = _mpi
    OBJ_DIR     = bin/mpi
    USE_CUDA_FLAG = 0
    USE_MPI_FLAG  = 1
else
    CXX       = g++
    CXXFLAGS  = -std=c++20 -Wall -O3 -I./include
    MODE_SUFFIX = _seq
    OBJ_DIR     = bin/seq
    USE_CUDA_FLAG = 0
    USE_MPI_FLAG  = 0
endif

# ---- Auto-detect CUDA architecture (fallback to sm_70) -------------------
ifeq ($(USE_CUDA_FLAG),1)
    CUDA_ARCH ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.')
    ifeq ($(CUDA_ARCH),)
        CUDA_ARCH = 70
    endif
    NVCCFLAGS += -arch=sm_$(CUDA_ARCH)
endif

# ---- Directories ----------------------------------------------------------
SRC_DIR = src
INC_DIR = include
BIN_DIR = bin

LDFLAGS_COMMON = -lm

ifeq ($(MODE),cuda_mpi)
    LDFLAGS_COMMON += $(MPI_LDFLAGS)
endif

# Add CUDA runtime link flags when needed
ifeq ($(USE_CUDA_FLAG),1)
    CUDA_LIB_DIR ?= $(shell dirname $$(which nvcc 2>/dev/null) 2>/dev/null)/../lib64
    LDFLAGS_COMMON += -L$(CUDA_LIB_DIR) -lcudart
endif

# ---- Core C++ sources (always compiled with CXX) -------------------------
CORE_CPP_SOURCES = \
    $(SRC_DIR)/DatasetLoader.cpp \
    $(SRC_DIR)/Simulation.cpp \
    $(SRC_DIR)/PerformanceLogger.cpp \
    $(SRC_DIR)/NaiveSimulation.cpp \
    $(SRC_DIR)/MpiNaiveSimulation.cpp \
    $(SRC_DIR)/BarnesHutSimulation.cpp \
    $(SRC_DIR)/CheckpointManager.cpp

CORE_CPP_OBJECTS = $(CORE_CPP_SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# ---- CUDA sources (only when USE_CUDA_FLAG=1) ----------------------------
ifeq ($(USE_CUDA_FLAG),1)
    CORE_CU_SOURCES = \
        $(SRC_DIR)/CudaNaiveSimulation.cu
    CORE_CU_OBJECTS = $(CORE_CU_SOURCES:$(SRC_DIR)/%.cu=$(OBJ_DIR)/%.cu.o)
else
    CORE_CU_SOURCES =
    CORE_CU_OBJECTS =
endif

CORE_OBJECTS = $(CORE_CPP_OBJECTS) $(CORE_CU_OBJECTS)

# ---- App sources ----------------------------------------------------------
HEADLESS_SRC = $(SRC_DIR)/apps/headless_main.cpp
HEADLESS_OBJ = $(HEADLESS_SRC:$(SRC_DIR)/apps/%.cpp=$(OBJ_DIR)/apps/%.o)

HEADLESS_TARGET = $(BIN_DIR)/nbody_headless$(MODE_SUFFIX)
VTK_CONVERTER   = $(BIN_DIR)/checkpoint_to_vtk

# ---- Targets --------------------------------------------------------------

all: $(BIN_DIR) $(HEADLESS_TARGET) $(VTK_CONVERTER)
	@echo "Built MODE=$(MODE) → $(HEADLESS_TARGET)"

$(BIN_DIR):
	mkdir -p $(BIN_DIR) $(OBJ_DIR) $(OBJ_DIR)/apps $(OBJ_DIR)/renderers

# Link headless binary
ifeq ($(USE_CUDA_FLAG),1)
$(HEADLESS_TARGET): $(CORE_OBJECTS) $(HEADLESS_OBJ)
	$(NVCC) $(NVCCFLAGS) -o $@ $^ $(LDFLAGS_COMMON)
	@echo "Linked $@ (CUDA enabled)"
else
$(HEADLESS_TARGET): $(CORE_OBJECTS) $(HEADLESS_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS_COMMON)
	@echo "Linked $@"
endif

# VTK converter (standalone, no deps on sim code)
$(VTK_CONVERTER): $(SRC_DIR)/tools/checkpoint_to_vtk.cpp
	@mkdir -p $(BIN_DIR)
	g++ -std=c++17 -Wall -O2 -I./include -o $@ $<
	@echo "Built $@"

# ---- Compile rules --------------------------------------------------------

# C++ → .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# CUDA → .cu.o
$(OBJ_DIR)/%.cu.o: $(SRC_DIR)/%.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

# ---- Phony targets --------------------------------------------------------

.PHONY: clean run debug help seq mpi cuda cuda_mpi vtk

clean:
	rm -rf $(BIN_DIR)

seq:
	$(MAKE) MODE=seq

mpi:
	$(MAKE) MODE=mpi

cuda:
	$(MAKE) MODE=cuda

cuda_mpi:
	$(MAKE) MODE=cuda_mpi

vtk: $(VTK_CONVERTER)

run: $(HEADLESS_TARGET)
	./$(HEADLESS_TARGET)

debug: CXXFLAGS += -g -DDEBUG
debug: clean all

help:
	@echo "N-Body Simulation — Build Targets"
	@echo ""
	@echo "  make seq          - Sequential CPU           (default)"
	@echo "  make mpi          - MPI distributed CPU"
	@echo "  make cuda         - Single-GPU CUDA"
	@echo "  make cuda_mpi     - Multi-GPU CUDA + MPI"
	@echo "  make vtk          - VTK converter only"
	@echo "  make clean        - Remove all build artifacts"
	@echo "  make run          - Build & run headless"
	@echo "  make debug        - Debug build"
	@echo ""
	@echo "Outputs:"
	@echo "  bin/nbody_headless_seq       - Sequential binary"
	echo "  bin/nbody_headless_mpi       - MPI binary"
	@echo "  bin/nbody_headless_cuda      - CUDA binary"
	@echo "  bin/nbody_headless_cuda_mpi  - CUDA+MPI binary"
	@echo "  bin/checkpoint_to_vtk        - Checkpoint → VTK"
	@echo "Override CUDA arch:  make cuda CUDA_ARCH=86  (for sm_86)"
	@echo ""
	@echo "VTK Converter Usage:"
	@echo "  ./bin/checkpoint_to_vtk simulation_output.bin ./vtk_output [--every N]"