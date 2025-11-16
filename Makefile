# Configuration: set to 'mpi' or 'seq'
MODE ?= seq

ifeq ($(MODE),mpi)
    CXX = mpic++
    CXXFLAGS = -std=c++20 -Wall -O0 -I./include -DUSE_MPI
    MODE_SUFFIX = _mpi
else
    CXX = g++
    CXXFLAGS = -std=c++20 -Wall -O0 -I./include
    MODE_SUFFIX = _seq
endif

UNAME_S := $(shell uname -s)

LDFLAGS_COMMON = -lm
LDFLAGS_GLUT   =

ifeq ($(UNAME_S),Linux)
    LDFLAGS_GLUT += -lglut -lGLU -lGL
endif

ifeq ($(UNAME_S),Darwin)
    LDFLAGS_GLUT += -framework OpenGL -framework GLUT
endif

SRC_DIR = src
INC_DIR = include
BIN_DIR = bin

# ----- core sources (no GL) -----
CORE_SOURCES = \
    $(SRC_DIR)/DatasetLoader.cpp \
    $(SRC_DIR)/Simulation.cpp \
    $(SRC_DIR)/PerformanceLogger.cpp \
    $(SRC_DIR)/NaiveSimulation.cpp \
    $(SRC_DIR)/MpiNaiveSimulation.cpp \
    $(SRC_DIR)/BarnesHutSimulation.cpp \
    $(SRC_DIR)/CheckpointManager.cpp

CORE_OBJECTS = $(CORE_SOURCES:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)

# ----- app sources -----
HEADLESS_SRC = $(SRC_DIR)/apps/headless_main.cpp
HEADLESS_OBJ = $(HEADLESS_SRC:$(SRC_DIR)/apps/%.cpp=$(BIN_DIR)/apps/%.o)

VIEWER_SRC = \
    $(SRC_DIR)/apps/viewer_main.cpp \
    $(SRC_DIR)/renderers/GlutRenderer.cpp

VIEWER_OBJS = $(VIEWER_SRC:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)

# targets with mode suffix
HEADLESS_TARGET = $(BIN_DIR)/nbody_headless$(MODE_SUFFIX)
VIEWER_TARGET   = $(BIN_DIR)/nbody_viewer$(MODE_SUFFIX)

all: $(BIN_DIR) $(HEADLESS_TARGET) $(VIEWER_TARGET)
	@echo "Built $(MODE) version"

$(BIN_DIR):
	mkdir -p $(BIN_DIR) $(BIN_DIR)/apps $(BIN_DIR)/renderers

# link
$(HEADLESS_TARGET): $(CORE_OBJECTS) $(HEADLESS_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS_COMMON)
	@echo "Built $@"

$(VIEWER_TARGET): $(CORE_OBJECTS) $(VIEWER_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS_COMMON) $(LDFLAGS_GLUT)
	@echo "Built $@"

# compile pattern rules
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: clean run run_viewer debug help seq mpi
clean:
	rm -rf $(BIN_DIR)

# Convenience targets
seq:
	$(MAKE) MODE=seq

mpi:
	$(MAKE) MODE=mpi

run: $(HEADLESS_TARGET)
	./$(HEADLESS_TARGET)

run_viewer: $(VIEWER_TARGET)
	./$(VIEWER_TARGET)

debug: CXXFLAGS += -g -DDEBUG
debug: clean all

help:
	@echo "Targets:"
	@echo "  make seq          - Build sequential version (default)"
	@echo "  make mpi          - Build MPI version"
	@echo "  make MODE=mpi     - Alternative way to build MPI"
	@echo "  all               - Build headless + viewer"
	@echo "  run               - Run headless app"
	@echo "  run_viewer        - Run OpenGL viewer"
	@echo "  debug             - Debug build"
	@echo "  clean             - Clean artifacts"
	@echo ""
	@echo "Examples:"
	@echo "  make seq          - builds bin/nbody_headless_seq"
	@echo "  make mpi          - builds bin/nbody_headless_mpi"