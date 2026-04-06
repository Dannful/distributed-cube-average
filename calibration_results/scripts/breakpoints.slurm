#!/bin/bash
#SBATCH --job-name=smpi-breakpoints
#SBATCH --partition=PARTITION_NAME
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE

OUTPUT_FILE="data/breakpoints"

BP1=$(mpirun -np $SLURM_NTASKS -machinefile $MACHINEFILE --mca btl ^openib nixw nix develop --command bp_search1 | tail -1)
BP2=$(mpirun -np $SLURM_NTASKS -machinefile $MACHINEFILE --mca btl ^openib nixw nix develop --command bp_search2 | tail -1)

echo "$BP1" > $OUTPUT_FILE
echo "$BP2" >> $OUTPUT_FILE

echo "Breakpoints written to $OUTPUT_FILE:"
cat $OUTPUT_FILE
