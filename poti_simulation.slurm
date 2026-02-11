#!/bin/bash
#SBATCH --job-name=distributed-cube-distributed
#SBATCH --partition=poti
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

nixw nix run .#run-poti-experiments
