#!/bin/bash

CUBE_SIZE_X=100
CUBE_SIZE_Y=100
CUBE_SIZE_Z=100
NUM_ITERATIONS=6
NUM_PROCESSES=16
STENCIL_SIZE=5

make all
mpirun --map-by :OVERSUBSCRIBE -np $NUM_PROCESSES ./bin/distributed-cube-average $CUBE_SIZE_X $CUBE_SIZE_Y $CUBE_SIZE_Z $NUM_ITERATIONS $STENCIL_SIZE ./predicted.dc
mpirun --map-by :OVERSUBSCRIBE -np 1 ./bin/distributed-cube-average $CUBE_SIZE_X $CUBE_SIZE_Y $CUBE_SIZE_Z $NUM_ITERATIONS $STENCIL_SIZE ./ground_truth.dc
Rscript CompareResults.R
