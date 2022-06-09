#!/bin/bash
#BSUB -J valgrind-problem2
#BSUB -q ser
#BSUB -n 1


module purge
module load intel/2018.4
module load mpi/intel/2018.4
module load valgrind/3.14.0

valgrind mpirun ./implicit -h 0.01 -dt 0.0000001 -maxit 100 -restart 0 > implicit-valgrind.log 2>&1
