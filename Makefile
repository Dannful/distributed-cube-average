CC = mpicc

SRCDIR = src
INCDIR = include
BUILDDIR = bin
OBJDIR = $(BUILDDIR)/obj

CFLAGS = -I$(INCDIR) -Wall -g -laky

SOURCES = $(wildcard $(SRCDIR)/*.c)

OBJECTS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES))

TARGET = $(BUILDDIR)/distributed-cube-average

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@echo "Cleaning compilation files..."
	rm -rf $(BUILDDIR)

.PHONY: all clean
