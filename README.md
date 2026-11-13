# Distributed cube average

This repository contains a distributed 3D seismic wave simulator that solves the Fletcher
equations. It is parallelized with MPI (3D Cartesian decomposition with halo exchange) and has
pluggable compute backends (OpenMP or CUDA). It can also be simulated at scale with
[SimGrid](https://simgrid.org/).

## Building

### Nix

In order to build using Nix, it is sufficient to run:
```bash
nix build
```
which will create a `result` symlink pointing to the build output in the nix store.

### Makefile

#### Dependencies

- An MPI implementation providing `mpicc`/`mpirun` (e.g. [OpenMPI](https://www.open-mpi.org/)
  or MPICH) and a C compiler with OpenMP support (e.g. GCC). Argument parsing uses `argp`,
  which is part of glibc. On Debian/Ubuntu:
  ```bash
  sudo apt install build-essential openmpi-bin libopenmpi-dev
  ```
- For the `cuda` and `cuda_aware` backends: the CUDA toolkit (`nvcc`), with `CUDA_HOME` set
  (e.g. `/usr/local/cuda`). `cuda_aware` also requires an MPI built with CUDA support.
- For the `simgrid` and `simgrid_cuda` backends: [SimGrid](https://simgrid.org/), which
  provides `smpicc`/`smpirun`.
- Optional profilers: [mpiP](https://github.com/LLNL/mpiP) (`libmpiP`) or
  [akypuera](https://github.com/schnorr/akypuera) (`libaky`).

#### Compiling

```bash
make                                   # OpenMP backend, no profiling
make BACKEND=openmp PROFILE=mpip       # with mpiP profiling
make BACKEND=cuda ARCH=sm_80           # CUDA, host-staged halos (default ARCH is sm_89)
make BACKEND=cuda_aware ARCH=sm_80     # CUDA, halos sent directly from GPU memory
make BACKEND=simgrid                   # OpenMP backend compiled for SimGrid
```

Supported values are `BACKEND=openmp|cuda|cuda_aware|simgrid|simgrid_cuda` and
`PROFILE=none|mpip|akypuera`. The executable is `bin/dc`.

Run `make clean` when switching backends: all backends share the same `bin/obj/` directory and
stale objects are not rebuilt.

## Running

### Program arguments

All arguments are mandatory:

| Option                         | Meaning                                                   |
|--------------------------------|-----------------------------------------------------------|
| `--size-x`, `--size-y`, `--size-z` | Physical domain size, in grid points                  |
| `--dx`, `--dy`, `--dz`         | Spatial step                                              |
| `--dt`                         | Time step                                                 |
| `--time-max`, `-t`             | Simulated time (iterations = `ceil(time-max / dt)`)       |
| `--absorption`, `-a`           | Width of the absorption zone, in grid points              |
| `--output-file`, `-o`          | Path of the binary output file                            |

Any number of MPI processes can be used: `MPI_Dims_create` factors them into a 3D topology, and
rank 0 acts as the coordinator.

### Nix

```bash
# nix run . -- BACKEND NP [program arguments]   (BACKEND: openmp | cuda | cuda_aware)
nix run . -- openmp 4 --size-x=100 --size-y=100 --size-z=100 \
  --dx=0.01 --dy=0.01 --dz=0.01 --dt=0.001 --time-max=0.5 \
  --absorption=10 --output-file=output.bin

# SimGrid: nix run .#simgrid -- NUM_HOSTS NET_BW NET_LAT [program arguments]
nix run .#simgrid -- 5 1Gbps 22.7us --size-x=100 ...
```

These wrappers also post-process traces and generate plots (see `nix/scripts.nix`).

### Manually (without Nix)

With real MPI (OpenMP backend):
```bash
export OMP_NUM_THREADS=4   # OpenMP threads per MPI rank
mpirun -np 4 --bind-to none ./bin/dc \
  --size-x 100 --size-y 100 --size-z 100 \
  --dx 0.01 --dy 0.01 --dz 0.01 \
  --dt 0.001 --time-max 0.5 \
  --absorption 10 \
  --output-file output.bin
```
`--bind-to none` (OpenMPI syntax) keeps OpenMPI from pinning each rank to a single core, which
would make all its OpenMP threads share that core.

The CUDA backends are run the same way (`mpirun -np N ./bin/dc ...`), typically one rank per
GPU.

With SimGrid (after `make BACKEND=simgrid`):
```bash
smpirun -np 5 -platform platform.xml -hostfile hostfile.txt ./bin/dc \
  --size-x 100 --size-y 100 --size-z 100 \
  --dx 0.01 --dy 0.01 --dz 0.01 \
  --dt 0.001 --time-max 0.5 \
  --absorption 10 \
  --output-file output.bin
```
The provided `platform.xml`/`hostfile.txt` describe 5 hosts (`host0`–`host4`); edit both to
simulate more processes. Add `-trace --cfg=tracing/filename:dc.trace` to obtain a Paje trace.
The calibrated SMPI network parameters used in the experiments are in `nix/scripts.nix` and can
be copied to the `smpirun` command line.

### Traces and plots (without Nix)

When built with `PROFILE=akypuera`, each rank writes a `rastro-*.rst` file. To convert and plot
them (requires akypuera, [PajeNG](https://github.com/schnorr/pajeng) and R with tidyverse):
```bash
aky_converter rastro-*.rst > dc.trace
pj_dump -z -l 9 dc.trace | grep ^State > dc.csv
Rscript plot.R dc.csv
```

## LSP

### Nix

There is a [flake](./flake.nix) that can be used to set up the development environment with all necessary dependencies, including `clangd`. To use it, you can run:
```bash
nix develop
```

### Manually

It is recommended to use [clangd](https://github.com/clangd/clangd) as the Language Server Protocol (LSP) for this project, since a [.clangd](./.clangd) file has already been created for this purpose.
