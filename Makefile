CC = mpicc
NVCC = nvcc
BACKEND ?= openmp
ARCH ?= sm_89

# Try to detect SimGrid path if not provided, defaults to /usr
SIMGRID_PATH ?= /usr

SRCDIR = src
INCDIR = include
BUILDDIR = bin
OBJDIR = $(BUILDDIR)/obj

CFLAGS = -I$(INCDIR) -Wall -O3
LDFLAGS = -lm

# Standard NVCC flags
CUDA_CFLAGS = -I$(INCDIR) -gencode arch=compute_$(subst sm_,,$(ARCH)),code=$(ARCH) -allow-unsupported-compiler

CUDA_LINK_LIBS = -L/usr/local/cuda/lib64 -lcudart

# Generic sources, common to all backends
SOURCES_C_COMMON = $(filter-out $(wildcard $(SRCDIR)/*_propagate.c) $(SRCDIR)/device_data.c $(SRCDIR)/derivatives.c $(SRCDIR)/sample.c, $(wildcard $(SRCDIR)/*.c))

# 1. OpenMP or Standard SimGrid (CPU-only)
ifeq ($(BACKEND), $(filter $(BACKEND), openmp simgrid))
    ifeq ($(BACKEND), simgrid)
        CC := smpicc
        CFLAGS += -D_SMPI_ADD_SEMANTICS 
    endif
    SOURCES_C := $(SOURCES_C_COMMON) $(SRCDIR)/openmp_propagate.c $(SRCDIR)/device_data.c
    SOURCES_CUDA :=
    CFLAGS += -fopenmp
    LDFLAGS += -fopenmp

# 2. Native CUDA (Default MPI + GPU)
else ifeq ($(BACKEND), cuda)
    SOURCES_C := $(SOURCES_C_COMMON)
    SOURCES_CUDA := $(SRCDIR)/cuda_propagate.cu $(SRCDIR)/device_data.cu
    LDFLAGS += $(CUDA_LINK_LIBS)

# 3. SimGrid + CUDA (Simulated MPI + GPU)
else ifeq ($(BACKEND), simgrid_cuda)
    CC := smpicc
    SOURCES_C := $(SOURCES_C_COMMON)
    SOURCES_CUDA := $(SRCDIR)/cuda_propagate.cu $(SRCDIR)/device_data.cu
    
    CUDA_CFLAGS += -I$(SIMGRID_PATH)/include
    CUDA_CFLAGS += -ccbin g++ 
    CFLAGS += -D_SMPI_ADD_SEMANTICS
    CUDA_CFLAGS += -D_SMPI_ADD_SEMANTICS
    
    LDFLAGS += $(CUDA_LINK_LIBS)

else
    $(error Unsupported backend: $(BACKEND). Supported: 'openmp', 'simgrid', 'cuda', 'simgrid_cuda')
endif

OBJECTS_C = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES_C))
OBJECTS_CUDA = $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(SOURCES_CUDA))
OBJECTS = $(OBJECTS_C) $(OBJECTS_CUDA)

TARGET = $(BUILDDIR)/dc

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJECTS_C): $(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(@D)
	$(CC) -c $< -o $@ $(CFLAGS)

$(OBJECTS_CUDA): $(OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) -c $< -o $@ $(CUDA_CFLAGS)

clean:
	@echo "Cleaning compilation files..."
	rm -rf $(BUILDDIR)

.PHONY: all clean
