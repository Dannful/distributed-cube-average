#!/bin/bash
#SBATCH --job-name=distributed-cube-distributed
#SBATCH --partition=draco
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

MACHINEFILE="nodes.$SLURM_JOB_ID"
srun -l hostname | sort -n | awk '{print $2}' > $MACHINEFILE

nix build .#default
mpirun -np $SLURM_NTASKS \
	-machinefile $MACHINEFILE \
	--mca btl ^openib \
	--mca btl_tcp_if_include eno2 \
	--bind-to none \
	nix shell -c ./result/bin/distributed-cube-average --size-x=100 --size-y=100 --size-z=100 --absorption=80 --dx=3 --dy=3 --dz=3 --dt=0.000110 --time-max=30 --output-file=./validation/predicted.dc
