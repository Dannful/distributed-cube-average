#!/bin/bash
#SBATCH --job-name=distributed-cube-distributed
#SBATCH --partition=cei
#SBATCH --nodes=6
#SBATCH --ntasks=6
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE

nix build .#dc

for i in 32 64 128 256 512; do
  i=$((i - 12))
  mpirun -np $SLURM_NTASKS \
    -machinefile $MACHINEFILE \
    --mca btl ^openib \
    --mca btl_tcp_if_include eno2 \
    --bind-to none \
    nix shell -c ./result/bin/dc openmp mpip --size-x=$i --size-y=$i --size-z=$i --absorption=2 --dx=1e-1 --dy=1e-1 --dz=1e-1 --dt=1e-6 --time-max=1e-4 --output-file=./validation/predicted.dc >execution.log 2>&1

  # Extract mpiP report from the log (everything from "@ mpiP" to the end)
  sed -n '/^@ mpiP/,$p' execution.log >dc.mpiP

  path="analysis/$(date +'%Y-%m-%d')/$((i + 12))"
  mkdir -p "$path"
  mv *.mpiP "$path" 2>/dev/null
  rm -f execution.log
done
