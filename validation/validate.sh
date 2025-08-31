#!/bin/bash

USE_SEQUENTIAL="true"
USE_DISTRIBUTED="true"

CUBE_SIZE_X=100
CUBE_SIZE_Y=100
CUBE_SIZE_Z=100
NUM_ITERATIONS=6
NUM_PROCESSES=8
STENCIL_SIZE=3

make all
if [[ "$USE_DISTRIBUTED" == "true" ]]; then
  mpirun --map-by :OVERSUBSCRIBE -np $NUM_PROCESSES ../bin/distributed-cube-average $CUBE_SIZE_X $CUBE_SIZE_Y $CUBE_SIZE_Z $NUM_ITERATIONS $STENCIL_SIZE ./predicted.dc
fi
if [[ "$USE_SEQUENTIAL" == "true" ]]; then
  mpirun --map-by :OVERSUBSCRIBE -np 1 ../bin/distributed-cube-average $CUBE_SIZE_X $CUBE_SIZE_Y $CUBE_SIZE_Z $NUM_ITERATIONS $STENCIL_SIZE ./ground_truth.dc
fi
if [[ "$USE_SEQUENTIAL" == "true" && "$USE_DISTRIBUTED" == "true" ]]; then
  Rscript CompareResults.R
fi
