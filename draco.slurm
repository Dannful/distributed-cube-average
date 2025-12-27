#!/bin/bash
#SBATCH --job-name=distributed-cube-distributed
#SBATCH --partition=draco
#SBATCH --nodes=6
#SBATCH --ntasks=6
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE
nix run .#draco -- $SLURM_NTASKS $MACHINEFILE
