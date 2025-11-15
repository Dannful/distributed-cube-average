CC = mpicc
NVCC = nvcc
BACKEND ?= openmp
ARCH ?= sm_89

SRCDIR = src
INCDIR = include
BUILDDIR = bin
OBJDIR = $(BUILDDIR)/obj

CFLAGS = -I$(INCDIR) -Wall -O3
LDFLAGS = -lm
CUDA_CFLAGS = -I$(INCDIR) -gencode arch=compute_$(subst sm_,,$(ARCH)),code=$(ARCH) -allow-unsupported-compiler
CUDA_LDFLAGS = -L/usr/local/cuda/lib64 -lcudart

# Generic sources, common to all backends
SOURCES_C_COMMON = $(filter-out $(wildcard $(SRCDIR)/*_propagate.c) $(SRCDIR)/device_data.c $(SRCDIR)/derivatives.c $(SRCDIR)/sample.c, $(wildcard $(SRCDIR)/*.c))

# Add backend-specific sources and flags
ifeq ($(BACKEND), openmp)
	SOURCES_C := $(SOURCES_C_COMMON) $(SRCDIR)/openmp_propagate.c $(SRCDIR)/device_data.c
	SOURCES_CUDA :=
	CFLAGS += -fopenmp
	LDFLAGS += -fopenmp
else ifeq ($(BACKEND), cuda)
	SOURCES_C := $(SOURCES_C_COMMON)
	SOURCES_CUDA := $(SRCDIR)/cuda_propagate.cu $(SRCDIR)/device_data.cu
else
	$(error Unsupported backend: $(BACKEND). Currently, only 'openmp' or 'cuda' are supported)
endif

OBJECTS_C = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES_C))
OBJECTS_CUDA = $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(SOURCES_CUDA))
OBJECTS = $(OBJECTS_C) $(OBJECTS_CUDA)

TARGET = $(BUILDDIR)/dc

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
ifeq ($(BACKEND), cuda)
	$(CC) -o $@ $^ $(LDFLAGS) $(CUDA_LDFLAGS)
else
	$(CC) -o $@ $^ $(LDFLAGS)
endif

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
