#!/bin/bash
#SBATCH --job-name=smpi-breakpoints
#SBATCH --partition=PARTITION_NAME
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE

OUTPUT_FILE="data/breakpoints"

BP_SEARCH1_CMD="nixw nix develop --command bp_search1"
BP_SEARCH2_CMD="nixw nix develop --command bp_search2"

nixw nix develop --command python3 scripts/find_breakpoints.py \
  "$BP_SEARCH1_CMD" \
  "$BP_SEARCH2_CMD" \
  -np $SLURM_NTASKS \
  -machinefile $MACHINEFILE \
  --mca btl ^openib \
  --mca btl_tcp_if_include 192.168.30.0/24 \
  > breakpoints_search.log 2>&1

if [ $? -ne 0 ]; then
  echo "ERROR: Breakpoint search failed. See breakpoints_search.log for details."
  cat breakpoints_search.log
  exit 1
fi

BP1=$(grep "First breakpoint:" breakpoints_search.log | awk '{print $NF}')
BP2=$(grep "Second breakpoint:" breakpoints_search.log | awk '{print $NF}')

echo "$BP1" > $OUTPUT_FILE
echo "$BP2" >> $OUTPUT_FILE

echo "Breakpoints written to $OUTPUT_FILE:"
cat $OUTPUT_FILE
