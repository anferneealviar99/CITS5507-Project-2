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

export OMP_NUM_THREADS=2
mpicc -fopenmp -o mpi_expr mpi_experiment.c
srun ./mpi_expr
