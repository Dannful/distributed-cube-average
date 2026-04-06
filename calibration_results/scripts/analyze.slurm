#!/bin/bash
#SBATCH --job-name=smpi-analyze
#SBATCH --partition=PARTITION_NAME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

mkdir -p outputs

nixw nix develop --command python src/mpi_calibration.py --data-dir ./data --output-dir ./outputs

echo "Analysis complete. Results in calibration_results/outputs/"
