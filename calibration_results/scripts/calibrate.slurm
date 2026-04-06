#!/bin/bash
#SBATCH --job-name=smpi-calibrate
#SBATCH --partition=PARTITION_NAME
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE

cd data

nixw nix develop ..#default --command python3 ../scripts/generate_experiment_plan.py exp_plan.csv

mpirun -np $SLURM_NTASKS \
  -machinefile ../$MACHINEFILE \
  --mca btl ^openib \
  nixw nix develop ..#default --command calibrate -d . -p exp -s exp_plan.csv

echo "Calibration complete. Data written to calibration_results/data/"
