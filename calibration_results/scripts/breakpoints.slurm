#!/bin/bash
#SBATCH --job-name=smpi-breakpoints
#SBATCH --partition=poti
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE

BP_SEARCH1_CMD="nixw nix develop --command bp_search1"
BP_SEARCH2_CMD="nixw nix develop --command bp_search2"

OUTPUT_FILE="../data/breakpoints"

nixw nix develop --command python3 scripts/find_breakpoints.py \
  "$BP_SEARCH1_CMD" \
  "$BP_SEARCH2_CMD" \
  -np $SLURM_NTASKS \
  -machinefile $MACHINEFILE \
  --mca pml ob1 \
  --mca btl tcp,self \
  --mca btl_tcp_if_include enp4s0 \
  --mca oob_tcl_if_include enp4s0 \
  > breakpoints_search.log 2>&1

if [ $? -ne 0 ]; then
  echo "ERROR: Breakpoint search failed. See breakpoints_search.log for details."
  cat breakpoints_search.log
  exit 1
fi

# Extract last word of each line; "no breakpoint found" lines produce a non-integer
BP1_RAW=$(grep "First breakpoint:" breakpoints_search.log | awk '{print $NF}')
BP2_RAW=$(grep "Second breakpoint:" breakpoints_search.log | awk '{print $NF}')

is_integer() { [[ "$1" =~ ^[0-9]+$ ]]; }

if is_integer "$BP2_RAW"; then
  BP_DEADLOCK="$BP2_RAW"
else
  echo "WARNING: No deadlock breakpoint found. Defaulting to 65536."
  BP_DEADLOCK=65536
fi

# Corruption threshold absent → use deadlock threshold for both
if is_integer "$BP1_RAW"; then
  BP_CORRUPTION="$BP1_RAW"
else
  echo "WARNING: No corruption breakpoint found. Using deadlock threshold ($BP_DEADLOCK) for both."
  BP_CORRUPTION="$BP_DEADLOCK"
fi

printf '%s\n%s\n' "$BP_CORRUPTION" "$BP_DEADLOCK" > "$OUTPUT_FILE"
echo "Breakpoints written to $OUTPUT_FILE:"
cat "$OUTPUT_FILE"
