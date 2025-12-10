# Configuration: set to 'mpi' or 'seq'
MODE ?= seq

ifeq ($(MODE),mpi)
    CXX = mpic++
    CXXFLAGS = -std=c++20 -Wall -O3 -I./include -DUSE_MPI
    MODE_SUFFIX = _mpi
    OBJ_DIR = bin/mpi
else
    CXX = g++
    CXXFLAGS = -std=c++20 -Wall -O3 -I./include
    MODE_SUFFIX = _seq
    OBJ_DIR = bin/seq
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

CORE_OBJECTS = $(CORE_SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# ----- app sources -----
HEADLESS_SRC = $(SRC_DIR)/apps/headless_main.cpp
HEADLESS_OBJ = $(HEADLESS_SRC:$(SRC_DIR)/apps/%.cpp=$(OBJ_DIR)/apps/%.o)

# Viewer always built without MPI
VIEWER_SRC = \
    $(SRC_DIR)/apps/viewer_main.cpp \
    $(SRC_DIR)/renderers/GlutRenderer.cpp

VIEWER_OBJS = $(VIEWER_SRC:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)

# Only headless gets mode suffix
HEADLESS_TARGET = $(BIN_DIR)/nbody_headless$(MODE_SUFFIX)
VIEWER_TARGET   = $(BIN_DIR)/nbody_viewer
VTK_CONVERTER   = $(BIN_DIR)/checkpoint_to_vtk

all: $(BIN_DIR) $(HEADLESS_TARGET) $(VIEWER_TARGET) $(VTK_CONVERTER)
	@echo "Built $(MODE) headless, viewer, and vtk converter"

$(BIN_DIR):
	mkdir -p $(BIN_DIR) $(OBJ_DIR) $(OBJ_DIR)/apps $(OBJ_DIR)/renderers $(BIN_DIR)/apps $(BIN_DIR)/renderers

# link headless (with mode-specific objects)
$(HEADLESS_TARGET): $(CORE_OBJECTS) $(HEADLESS_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS_COMMON)
	@echo "Built $@"

# link viewer (always sequential, recompile core objects without MPI)
$(VIEWER_TARGET): $(BIN_DIR)/seq_core.stamp $(VIEWER_OBJS)
	g++ -std=c++20 -Wall -O3 -I./include -o $@ \
		$(BIN_DIR)/viewer_objs/*.o $(VIEWER_OBJS) \
		$(LDFLAGS_COMMON) $(LDFLAGS_GLUT)
	@echo "Built $@"

# VTK converter (standalone, no dependencies on simulation code)
$(VTK_CONVERTER): $(SRC_DIR)/tools/checkpoint_to_vtk.cpp
	@mkdir -p $(BIN_DIR)
	g++ -std=c++17 -Wall -O2 -I./include -o $@ $<
	@echo "Built $@"

# Special rule to build core objects for viewer (always sequential)
$(BIN_DIR)/seq_core.stamp: $(CORE_SOURCES)
	@mkdir -p $(BIN_DIR)/viewer_objs
	@for src in $(CORE_SOURCES); do \
		obj=$$(echo $$src | sed 's|$(SRC_DIR)/|$(BIN_DIR)/viewer_objs/|; s|\.cpp$$|.o|'); \
		mkdir -p $$(dirname $$obj); \
		g++ -std=c++20 -Wall -O3 -I./include -c -o $$obj $$src; \
	done
	@touch $@

# compile pattern rules for mode-specific objects
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# compile viewer objects (always without MPI)
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	g++ -std=c++20 -Wall -O3 -I./include -c -o $@ $<

.PHONY: clean run run_viewer debug help seq mpi vtk

clean:
	rm -rf $(BIN_DIR)

# Convenience targets
seq:
	$(MAKE) MODE=seq

mpi:
	$(MAKE) MODE=mpi

vtk: $(VTK_CONVERTER)

run: $(HEADLESS_TARGET)
	./$(HEADLESS_TARGET)

run_viewer: $(VIEWER_TARGET)
	./$(VIEWER_TARGET)

debug: CXXFLAGS += -g -DDEBUG
debug: clean all

help:
	@echo "Targets:"
	@echo "  make seq          - Build sequential headless + viewer + vtk converter"
	@echo "  make mpi          - Build MPI headless + viewer + vtk converter"
	@echo "  make vtk          - Build only the VTK converter"
	@echo "  all               - Build headless + viewer + vtk converter"
	@echo "  run               - Run headless app"
	@echo "  run_viewer        - Run OpenGL viewer"
	@echo "  debug             - Debug build"
	@echo "  clean             - Clean artifacts"
	@echo ""
	@echo "Outputs:"
	@echo "  bin/nbody_headless_seq  - Sequential headless"
	@echo "  bin/nbody_headless_mpi  - MPI headless"
	@echo "  bin/nbody_viewer        - Viewer (always sequential)"
	@echo "  bin/checkpoint_to_vtk   - Checkpoint to VTK converter"
	@echo ""
	@echo "VTK Converter Usage:"
	@echo "  ./bin/checkpoint_to_vtk simulation_output.bin ./vtk_output [--every N]"