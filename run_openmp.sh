#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # number of processes = 1 
#SBATCH --cpus-per-task=40      # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

./openmp 100 100 4
./openmp 100 100 8
./openmp 100 100 12
./openmp 100 100 16
./openmp 100 100 20
./openmp 100 100 24
./openmp 100 100 28
./openmp 100 100 32
/.openmp 100 100 36

#./openmp 1000 100 4
#./openmp 1000 100 8
#./openmp 1000 100 12
#./openmp 1000 100 16
#./openmp 1000 100 20
#./openmp 1000 100 24
#./openmp 1000 100 28
#./openmp 1000 100 32
#./openmp 1000 100 36

#./openmp 10000 100 4
#./openmp 10000 100 8
#./openmp 10000 100 12
#./openmp 10000 100 16
#./openmp 10000 100 20
#./openmp 10000 100 24
#./openmp 10000 100 28
#./openmp 10000 100 32
#./openmp 10000 100 36
