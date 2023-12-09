CC := g++-10
NVCC := nvcc
CFLAGS := -std=c++17 -g -fpic -MMD -Wall -fopenmp
NVCC_FLAGS := -std=c++17 -fpic -MMD -Xcudafe="--diag_suppress=20012 --diag_suppress=20236" -Xcompiler -fopenmp

INCLUDE := libraries/include
LIBS_DIR := libraries/lib
LIBS_DIR_GPU := /usr/local/cuda/lib64
LIBS :=
LIBS_GPU := cuda cudart cublas

TARGET = curv3d
SRC_DIR = src
BUILD_DIR = build

# Should not need to modify below.

CPU_BUILD_DIR = $(BUILD_DIR)/cpu
GPU_BUILD_DIR = $(BUILD_DIR)/gpu
PY_CPU_BUILD_DIR = $(BUILD_DIR)/py_cpu

SRC = $(wildcard $(SRC_DIR)/*/*.cpp) $(wildcard $(SRC_DIR)/*.cpp)

# Get source files and object files.
GCC_SRC = $(filter-out %.cu.cpp ,$(SRC))
NVCC_SRC = $(filter %.cu.cpp, $(SRC))
GCC_OBJ = $(GCC_SRC:$(SRC_DIR)/%.cpp=%.o)
NVCC_OBJ = $(NVCC_SRC:$(SRC_DIR)/%.cpp=%.o)

# If compiling for CPU, all go to GCC. Otherwise, they are split.
CPU_OBJ = $(addprefix $(CPU_BUILD_DIR)/,$(GCC_OBJ)) $(addprefix $(CPU_BUILD_DIR)/,$(NVCC_OBJ))
GPU_GCC_OBJ = $(addprefix $(GPU_BUILD_DIR)/,$(GCC_OBJ))
GPU_NVCC_OBJ = $(addprefix $(GPU_BUILD_DIR)/,$(NVCC_OBJ))

PY_CPU_OBJ = $(addprefix $(PY_CPU_BUILD_DIR)/,$(GCC_OBJ)) $(addprefix $(PY_CPU_BUILD_DIR)/,$(NVCC_OBJ))

# $(info $$GCC_SRC is [${GCC_SRC}])
# $(info $$NVCC_SRC is [${NVCC_SRC}])
# $(info $$GCC_OBJ is [${GCC_OBJ}])
# $(info $$NVCC_OBJ is [${NVCC_OBJ}])

# $(info $$CPU_OBJ is [${CPU_OBJ}])
# $(info $$GPU_GCC_OBJ is [${GPU_GCC_OBJ}])
# $(info $$GPU_NVCC_OBJ is [${GPU_NVCC_OBJ}])

HEADER = $(wildcard $(SRC_DIR)/*/*.h) $(wildcard $(SRC_DIR)/*.h)
CPU_DEPS = $(wildcard $(CPU_BUILD_DIR)/*.d)
GPU_DEPS = $(wildcard $(GPU_BUILD_DIR)/*.d)

SPACE := $(subst ,, )
COMMA := ,

INC := $(INCLUDE:%=-I%)
LIB := $(LIBS_DIR:%=-L%)
RPATH := $(subst $(SPACE),$(COMMA),$(LIBS_DIR))
LIB_GPU := $(LIBS_DIR_GPU:%=-L%)
LD := $(LIBS:%=-l%)
LD_GPU := $(LIBS_GPU:%=-l%)

PY_EXT := $(shell python3-config --extension-suffix)
PY_INC := $(shell python3-config --includes)


# Reminder:
# $< = first prerequisite
# $@ = the target which matched the rule
# $^ = all prerequisites

.PHONY: all docs clean

all : cpu gpu python docs

cpu: $(TARGET)CPU
gpu: $(TARGET)GPU
python: $(TARGET)PY

docs: SHELL:=/bin/bash
docs:
	if [ ! -d "docs/.venv" ]; then python3 -m venv "docs/.venv"; fi
	( \
		source docs/.venv/bin/activate;\
		pip install -r docs/requirements.txt; \
		doxygen docs/Doxyfile; \
		sphinx-build -b html docs/source docs/build/html; \
		deactivate \
	)

# CPU ONLY
$(TARGET)CPU: $(CPU_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB) $(LD) -Wl,-rpath,$(RPATH)

$(CPU_BUILD_DIR)/%.o $(CPU_BUILD_DIR)/%.cu.o: $(SRC_DIR)/%.cpp | $(CPU_BUILD_DIR)
	$(CC) $(CFLAGS) -c -o $@ $< $(INC) $(PY_INC)

# GPU ONLY
# For GPU, we need to build the NVCC objects, the NVCC linked object, and the
# regular ones. Then, we link them all together.
$(TARGET)GPU: $(GPU_BUILD_DIR)/link.o $(GPU_GCC_OBJ) | $(GPU_BUILD_DIR)
	$(CC) -g -DCUDA $(CFLAGS) $(GPU_NVCC_OBJ) $^ -o $@ $(INC) $(LIB) $(LIB_GPU) $(LD) $(LD_GPU)  -Wl,-rpath,$(RPATH)

$(GPU_BUILD_DIR)/link.o: $(GPU_NVCC_OBJ) | $(GPU_BUILD_DIR)
	$(NVCC) --device-link $^ -o $@

$(GPU_BUILD_DIR)/%.cu.o: $(SRC_DIR)/%.cu.cpp | $(GPU_BUILD_DIR)
	$(NVCC) $(NVCC_FLAGS) -DCUDA -x cu --device-c -o $@ $< $(INC)

$(GPU_BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(GPU_BUILD_DIR)
	$(CC) $(CFLAGS) -g -DCUDA -c -o $@ $< $(INC)

# PYTHON CPU

$(TARGET)PY: $(PY_CPU_OBJ)
	$(CC) $(CFLAGS) $^ -shared -o $(TARGET)$(PY_EXT) $(LIB) -Wl,--whole-archive $(LD) -Wl,--no-whole-archive

$(PY_CPU_BUILD_DIR)/%.o $(PY_CPU_BUILD_DIR)/%.cu.o: $(SRC_DIR)/%.cpp | $(PY_CPU_BUILD_DIR)
	$(CC) $(CFLAGS) -DPYTHON -c -o $@ $< $(INC) $(PY_INC)

-include $(CPU_DEPS)
-include $(GPU_DEPS)

$(CPU_BUILD_DIR) $(GPU_BUILD_DIR) $(PY_CPU_BUILD_DIR):
	mkdir -p $@

clean:
	rm -Rf $(BUILD_DIR) $(TARGET)* docs/build docs/.venv
