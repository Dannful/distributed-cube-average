#!/bin/bash
#SBATCH --job-name=distributed-cube-distributed
#SBATCH --partition=poti
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2 " slots=1"}' >$MACHINEFILE

nixw nix build .#dc

current_date=$(date +'%Y-%m-%d')

for i in {1..60}; do
  for j in 32 64 128 256 512; do
    j=$((j - 12))
    mpirun -np $SLURM_NTASKS \
      -machinefile $MACHINEFILE \
      --mca btl ^openib \
      --mca btl_tcp_if_include 192.168.30.0/24 \
      --bind-to none \
      nixw ./result/bin/dc cuda akypuera --size-x=$j --size-y=$j --size-z=$j --absorption=2 --dx=1e-1 --dy=1e-1 --dz=1e-1 --dt=1e-6 --time-max=1e-3 --output-file=./validation/predicted.dc
    path="analysis/$current_date/$i/$((j + 12))"
    mkdir -p "$path"
  done
done
