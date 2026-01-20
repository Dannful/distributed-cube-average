#!/bin/bash
#SBATCH --job-name=distributed-cube-distributed
#SBATCH --partition=poti
#SBATCH --nodes=5
#SBATCH --ntasks=5
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE

nix build .#dc

current_date=$(date +'%Y-%m-%d')

for i in {1..60}; do
  for j in 32 64 128 256 512; do
    j=$((j - 12))
    mpirun -np $SLURM_NTASKS \
      -machinefile $MACHINEFILE \
      --mca btl ^openib \
      --mca btl_tcp_if_include 192.168.30.0/24 \
      --bind-to none \
      nix shell -c ./result/bin/dc cuda mpip --size-x=$j --size-y=$j --size-z=$j --absorption=2 --dx=1e-1 --dy=1e-1 --dz=1e-1 --dt=1e-6 --time-max=1e-4 --output-file=./validation/predicted.dc >execution.log 2>&1

    # Extract mpiP report from the log (everything from "@ mpiP" to the end)
    sed -n '/^@ mpiP/,$p' execution.log >dc.mpiP

    path="analysis/$current_date/$i/$((j + 12))"
    mkdir -p "$path"
    mv *.mpiP "$path" 2>/dev/null
    rm -f execution.log
  done
done
