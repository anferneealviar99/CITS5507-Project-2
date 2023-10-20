#!/bin/sh
#SBATCH --account=courses0101
#SBATCH --partition=debug
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=4
#SBATCH --time=00:30:00
#SBATCH --exclusive
#SBATCH --mem-per-cpu=32G
#module load openmpi/4.0.5

export OMP_NUM_THREADS=16
mpicc -fopenmp -o fish_school_parallel fish_school_parallel.c
srun ./fish_school_parallel
