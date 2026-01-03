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
nix build .#cuda
for i in 32 64 128; do
  i=$((i - 12))
  mpirun -np $SLURM_NTASKS \
    -machinefile $MACHINEFILE \
    --mca btl ^openib \
    --mca btl_tcp_if_include 192.168.30.0/24 \
    --bind-to none \
    nix shell -c ./result/bin/dc --size-x=$i --size-y=$i --size-z=$i --absorption=2 --dx=1e-1 --dy=1e-1 --dz=1e-1 --dt=1e-6 --time-max=1e-4 --output-file=./validation/predicted.dc
  path="analysis/$(date +'%Y-%m-%d')/$((i + 12))"
  mkdir -p "$path"
  mv *.rst "$path"
done
