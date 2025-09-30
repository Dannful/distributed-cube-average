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
	nix shell -c ./result/bin/distributed-cube-average --size-x=1000 --size-y=1000 --size-z=1000 --absorption=800 --dx=3 --dy=3 --dz=3 --dt=0.000110 --time-max=30 --output-file=./validation/ground_truth.dc
