#!/bin/bash
#SBATCH --job-name=distributed-cube-sequential
#SBATCH --partition=draco
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

nix build .#default
mpirun -np $SLURM_NTASKS \
	nix shell -c ./result/bin/distributed-cube-average 1000 1000 1000 15 5 ./validation/ground_truth.dc
