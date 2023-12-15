#!/bin/bash

#SBATCH --nodes=1	     # number of nodes on which to run
#SBATCH --ntasks=1	     # number of tasks to run
#SBATCH --cpus-per-task=16    # number of cpus required per task
#SBATCH --partition=rome	     # partition requested
#SBATCH --time=48:00:00	     # wall-clock time limit

module load 2022
module load GCC/11.3.0
module load parallel/20220722-GCCcore-11.3.0

export OMP_NUM_THREADS=4
export GFORTRAN_UNBUFFERED_PRECONNECTED=y

file=d200mm_sensitivity_starikovskiy_0.03H2O_commands.txt

awk -v n=$SLURM_ARRAY_TASK_ID 'NR % 20 == n' $file | parallel --jobs 4 --ungroup
