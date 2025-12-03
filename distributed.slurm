#!/bin/bash
#SBATCH --job-name=distributed-cube-distributed
#SBATCH --partition=draco
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.slurm"
nix build .#cuda
mpirun -np $SLURM_NTASKS \
  -machinefile $MACHINEFILE \
  --mca btl ^openib \
  --mca btl_tcp_if_include eno2 \
  --bind-to none \
  nix shell -c ./result/bin/dc --size-x=500 --size-y=500 --size-z=500 --absorption=80 --dx=0.1 --dy=0.1 --dz=0.1 --dt=0.000110 --time-max=1 --output-file=./validation/predicted.dc
