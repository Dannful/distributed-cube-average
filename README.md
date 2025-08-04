# Distributed cube average

This repository contains an implementation of a simple average of a distributed cube using OpenMPI and C as a proof of concept for distributed computing.

## Building

### Nix

In order to build using Nix, it is sufficient to run:
```bash
nix build
```
which will create a `result` symlink pointing to the build output in the nix store.

### Makefile

First, ensure you have [OpenMPI](https://www.open-mpi.org/) installed on your system. You can then build the project using the following command:
```bash
make
```
which will output an executable named `distributed_cube_average` under the `bin/` directory.

## Running

### Nix

To run the program using Nix, you can execute the following command:
```bash
nix run
```

### Manually

To run the program, you can use the `mpirun` command. For example, to run the program with 4 processes, you can use:
```bash
mpirun -np 4 ./bin/distributed_cube_average
```

## LSP

### Nix

There is a [flake](./flake.nix) that can be used to set up the development environment with all necessary dependencies, including `clangd`. To use it, you can run:
```bash
nix develop
```

### Manually

It is recommended to use [clangd](https://github.com/clangd/clangd) as the Language Server Protocol (LSP) for this project, since a [.clangd](./.clangd) file has already been created for this purpose.
