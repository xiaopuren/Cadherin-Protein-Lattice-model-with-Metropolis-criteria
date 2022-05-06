#!/bin/sh
#$ -V
#$ -cwd
##########$ -q gpu.q  
#$ -q all.q@n411
########$ -pe mpi 1      ## Parallel jobs beyond 28 cores will require using mpi and the mpi parallel environment (-pe mpi X)
########$ -pe gpu-smp X  ## when X <= 20 if use multiple cores on one node
#$ -N ad

module load gcc

###/gs/gsfs0/users/zsu/KMC/TNF-3D/test/out
./out.o
