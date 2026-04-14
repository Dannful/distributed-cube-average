#!/bin/bash
#SBATCH --job-name=smpi-calibrate
#SBATCH --partition=poti
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE

cd data

echo "Prewarming..."
srun --ntasks=$SLURM_NTASKS --nodes=$SLURM_NNODES nixw nix build ..#platform-calibration --no-link
echo "Done."

echo "Generating experiments..."
nixw nix develop ..#default --command python3 ../scripts/generate_experiment_plan.py exp_plan.csv \
  --num-sizes 150 --repetitions 15 --max-size 134217728
echo "Done."

echo "Calibrating..."
mpirun -np $SLURM_NTASKS \
  -machinefile ../$MACHINEFILE \
  --mca btl ^openib \
  --mca btl_tcp_if_include 192.168.30.0/24 \
  --bind-to none \
  nixw nix shell ..#platform-calibration --command calibrate -m 0 -M 1000000000 -d . -p exp -s exp_plan.csv

echo "Calibration complete. Data written to calibration_results/data/"
