#!/bin/bash

#SBATCH --nodes=1	     # number of nodes on which to run
#SBATCH --ntasks=1	     # number of tasks to run
#SBATCH --cpus-per-task=128    # number of cpus required per task
#SBATCH --partition=rome	     # partition requested
#SBATCH --time=5-00:00:00	     # wall-clock time limit

module load 2022
module load GCC/11.3.0
module load parallel/20220722-GCCcore-11.3.0

export OMP_NUM_THREADS=1
export GFORTRAN_UNBUFFERED_PRECONNECTED=y

parallel --jobs 45 --ungroup < d200mm_streamer_comparison_commands.txt