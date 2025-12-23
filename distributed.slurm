#!/bin/bash
#SBATCH --job-name=distributed-cube-distributed
#SBATCH --partition=tupi
#SBATCH --nodes=6
#SBATCH --ntasks=6
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
srun -l hostname | sort -n | awk '{print $2}' >$MACHINEFILE
nix build .#cuda
mpirun -np $SLURM_NTASKS \
  -machinefile $MACHINEFILE \
  --mca btl ^openib \
  --mca btl_tcp_if_include 192.168.30.0/24 \
  --bind-to none \
  nix shell -c ./result/bin/dc --size-x=20 --size-y=20 --size-z=20 --absorption=2 --dx=0.1 --dy=0.1 --dz=0.1 --dt=1e-6 --time-max=1e-4 --output-file=./validation/predicted.dc
