CC = mpicc
BACKEND ?= openmp

SRCDIR = src
INCDIR = include
BUILDDIR = bin
OBJDIR = $(BUILDDIR)/obj

CFLAGS = -I$(INCDIR) -Wall -lm -laky

# Generic sources, excluding backend-specific implementations
SOURCES = $(filter-out $(wildcard $(SRCDIR)/*_propagate.c), $(wildcard $(SRCDIR)/*.c))

# Add backend-specific sources and flags
ifeq ($(BACKEND), openmp)
SOURCES += $(SRCDIR)/openmp_propagate.c
CFLAGS += -fopenmp
else ifeq ($(BACKEND), cuda)
$(error CUDA backend is not supported yet)
else
$(error Unsupported backend: $(BACKEND). Currently, only 'openmp' is supported)
endif

OBJECTS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES))

TARGET = $(BUILDDIR)/distributed-cube-average

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CC) -o $@ $^ $(CFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(@D)
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	@echo "Cleaning compilation files..."
	rm -rf $(BUILDDIR)

.PHONY: all clean
