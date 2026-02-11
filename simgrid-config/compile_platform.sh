#!/usr/bin/env bash
set -e

# Compile the S4U platform code in two modes:
# 1. As a shared library (.so) for smpirun to load.
# 2. As an executable to generate the hostfile.

if ! command -v pkg-config &> /dev/null; then
    echo "Error: pkg-config not found. Please enter the nix dev shell."
    exit 1
fi

SG_CFLAGS=$(pkg-config --cflags simgrid)
SG_LIBS=$(pkg-config --libs simgrid)

echo "SimGrid CFLAGS: $SG_CFLAGS"
echo "SimGrid LIBS:   $SG_LIBS"

echo "1. Compiling Shared Object (libplatform.so)..."
# We specifically use -std=c++17 as SimGrid often requires it
g++ -std=c++17 -shared -fPIC -o simgrid-config/libplatform.so simgrid-config/platform_s4u.cpp $SG_CFLAGS $SG_LIBS

echo "2. Compiling Generator Executable (generate_artifacts)..."
g++ -std=c++17 -o simgrid-config/generate_artifacts simgrid-config/platform_s4u.cpp $SG_CFLAGS $SG_LIBS

echo "Done."
echo ""
echo "================================================================================"
echo "                           USAGE GUIDE"
echo "================================================================================"
echo "Configuration is controlled via Environment Variables."
echo ""
echo "Required Variables:"
echo "  PLATFORM_NUM_HOSTS             # Number of hosts to simulate"
echo "  PLATFORM_NET_BW                # Network bandwidth (Host <-> Switch) (e.g. 1.25GBps)"
echo "  PLATFORM_NET_LAT               # Network latency (e.g. 50us)"
echo "  PLATFORM_GPU_BW                # PCIe/NVLink bandwidth (Host -> GPU) (e.g. 16GBps)"
echo "  PLATFORM_GPU_LAT               # PCIe/NVLink latency (e.g. 0us)"
echo "  PLATFORM_HOSTFILE              # Path for generated hostfile (e.g. simgrid-config/hostfile.txt)"
echo ""
echo "Step 1: Generate the hostfile"
echo "  export PLATFORM_NUM_HOSTS=8"
echo "  export PLATFORM_NET_BW=1.25GBps"
echo "  export PLATFORM_NET_LAT=50us"
echo "  export PLATFORM_GPU_BW=16GBps"
echo "  export PLATFORM_GPU_LAT=0us"
echo "  export PLATFORM_HOSTFILE=simgrid-config/hostfile.txt"
echo "  ./simgrid-config/generate_artifacts"
echo ""
echo "Step 2: Run the simulation"
echo "  smpirun \\"
echo "    -platform simgrid-config/libplatform.so \\"
echo "    -hostfile simgrid-config/hostfile.txt \\"
echo "    -np 8 \\"
echo "    ./your_mpi_application"
echo "================================================================================"