CC       = mpicc
NVCC     = nvcc
BACKEND ?= openmp
ARCH    ?= sm_89

SRCDIR   = src
INCDIR   = include
BUILDDIR = bin
OBJDIR   = $(BUILDDIR)/obj

CFLAGS   = -I$(INCDIR) -Wall -O3
LDFLAGS  = -lm

CUDA_CFLAGS    = -I$(INCDIR) -gencode arch=compute_$(subst sm_,,$(ARCH)),code=$(ARCH) -allow-unsupported-compiler
CUDA_LINK_LIBS = -L/usr/local/cuda/lib64 -lcudart

SOURCES_C_COMMON = $(filter-out $(wildcard $(SRCDIR)/*_propagate.c) $(SRCDIR)/device_data.c $(SRCDIR)/derivatives.c $(SRCDIR)/sample.c, $(wildcard $(SRCDIR)/*.c))

ifeq ($(BACKEND), openmp)
    SOURCES_C    := $(SOURCES_C_COMMON) $(SRCDIR)/openmp_propagate.c $(SRCDIR)/device_data.c
    SOURCES_CUDA :=
    CFLAGS       += -fopenmp
    LDFLAGS      += -fopenmp -laky

else ifeq ($(BACKEND), simgrid)
    CC           := smpicc
    SOURCES_C    := $(SOURCES_C_COMMON) $(SRCDIR)/openmp_propagate.c $(SRCDIR)/device_data.c
    SOURCES_CUDA :=
    CFLAGS       += -D_SMPI_ADD_SEMANTICS

else ifeq ($(BACKEND), cuda)
    SOURCES_C    := $(SOURCES_C_COMMON)
    SOURCES_CUDA := $(SRCDIR)/cuda_propagate.cu $(SRCDIR)/device_data.cu
    LDFLAGS      += $(CUDA_LINK_LIBS) -laky

else ifeq ($(BACKEND), simgrid_cuda)
    CC           := smpicc
    SOURCES_C    := $(SOURCES_C_COMMON)
    SOURCES_CUDA := $(SRCDIR)/cuda_propagate.cu $(SRCDIR)/device_data.cu
    
    CUDA_CFLAGS  += -Xcompiler -fPIC
    CFLAGS       += -D_SMPI_ADD_SEMANTICS
    CUDA_CFLAGS  += -D_SMPI_ADD_SEMANTICS
    CUDA_CFLAGS  += -ccbin g++
    LDFLAGS      += -L/usr/local/cuda/lib64 -lcudart_static -lstdc++ -lpthread -ldl -lrt

else
    $(error Unsupported backend: $(BACKEND). Supported: openmp, simgrid, cuda, simgrid_cuda)
endif

OBJECTS_C    = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES_C))
OBJECTS_CUDA = $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(SOURCES_CUDA))

TARGET = $(BUILDDIR)/dc

all: $(TARGET)

$(TARGET): $(OBJECTS_C) $(OBJECTS_CUDA)
	@mkdir -p $(@D)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJECTS_C): $(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(@D)
	$(CC) -c $< -o $@ $(CFLAGS)

$(OBJECTS_CUDA): $(OBJDIR)/%.o: $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) -c $< -o $@ $(CUDA_CFLAGS)

clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean
