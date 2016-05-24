#!/bin/sh
#PBS -N EIK3D
#PBS -l select=1:ncpus=12:ompthreads=12:mem=24gb,place=scatter:exclhost
#PBS -l walltime=06:00:00
#PBS -m n
#PBS -q bl2x220g7q

#cd $PBS_O_WORKDIR
echo "Working directory: $PBS_O_WORKDIR"
echo "Assigned nodes:"
#cat $PBS_NODEFILE
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo

./test.sh
