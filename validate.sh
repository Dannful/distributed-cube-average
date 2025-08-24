#!/bin/bash

CUBE_SIZE_X=4
CUBE_SIZE_Y=4
CUBE_SIZE_Z=4
NUM_ITERATIONS=1
NUM_PROCESSES=8
STENCIL_SIZE=1

make all
mpirun --map-by :OVERSUBSCRIBE -np $NUM_PROCESSES ./bin/distributed-cube-average $CUBE_SIZE_X $CUBE_SIZE_Y $CUBE_SIZE_Z $NUM_ITERATIONS $STENCIL_SIZE ./predicted.dc > 2
mpirun --map-by :OVERSUBSCRIBE -np 1 ./bin/distributed-cube-average $CUBE_SIZE_X $CUBE_SIZE_Y $CUBE_SIZE_Z $NUM_ITERATIONS $STENCIL_SIZE ./ground_truth.dc > 1
Rscript CompareResults.R
