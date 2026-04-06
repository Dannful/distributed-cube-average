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

nixw nix develop ..#default --command python3 -c "
import sys
n = int(sys.argv[1])
for i in range(1, n+1):
    print(int(1.1**i))
" 220 > sizes.txt

mpirun -np $SLURM_NTASKS \
  -machinefile ../$MACHINEFILE \
  --mca btl ^openib \
  nixw nix develop ..#default --command calibrate -d . -p exp -s sizes.txt

echo "Calibration complete. Data written to calibration_results/data/"
