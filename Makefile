CC = mpicc
NVCC = nvcc
BACKEND ?= openmp

SRCDIR = src
INCDIR = include
BUILDDIR = bin
OBJDIR = $(BUILDDIR)/obj

CFLAGS = -I$(INCDIR) -Wall -lm -laky
LDFLAGS = -lm -laky
CUDA_CFLAGS = -I$(INCDIR) --compiler-options "-Wall" -fmad=false --ftz=false --prec-div=true
CUDA_LDFLAGS = -L/usr/local/cuda/lib64 -lcudart

# Generic sources, excluding backend-specific implementations
SOURCES_C = $(filter-out $(wildcard $(SRCDIR)/*_propagate.c) $(SRCDIR)/derivatives.c $(SRCDIR)/sample.c, $(wildcard $(SRCDIR)/*.c))
SOURCES_CUDA =

# Add backend-specific sources and flags
ifeq ($(BACKEND), openmp)
	SOURCES_C += $(SRCDIR)/openmp_propagate.c
	CFLAGS += -fopenmp
	LDFLAGS += -fopenmp
else ifeq ($(BACKEND), cuda)
	SOURCES_CUDA += $(SRCDIR)/cuda_propagate.cu
else
	$(error Unsupported backend: $(BACKEND). Currently, only 'openmp' or 'cuda' are supported)
endif

OBJECTS_C = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES_C))
OBJECTS_CUDA = $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(SOURCES_CUDA))
OBJECTS = $(OBJECTS_C) $(OBJECTS_CUDA)

TARGET = $(BUILDDIR)/distributed-cube-average

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
ifeq ($(BACKEND), cuda)
	$(CC) -o $@ $^ $(LDFLAGS) $(CUDA_LDFLAGS)
else
	$(CC) -o $@ $^ $(LDFLAGS)
endif

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(@D)
	$(CC) -c $< -o $@ $(CFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) -c $< -o $@ $(CUDA_CFLAGS)

clean:
	@echo "Cleaning compilation files..."
	rm -rf $(BUILDDIR)

.PHONY: all clean
