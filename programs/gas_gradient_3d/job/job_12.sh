#!/bin/bash

#SBATCH --nodes=1	     # number of nodes on which to run
#SBATCH --ntasks=4	     # number of tasks to run
#SBATCH --cpus-per-task=32    # number of cpus required per task
#SBATCH --partition=thin	     # partition requested
#SBATCH --mem=120G             # memory requested
#SBATCH --time=3-00:00:00	     # wall-clock time limit


module load 2022
module load GCC/11.3.0

export OMP_NUM_THREADS=32
export GFORTRAN_UNBUFFERED_PRECONNECTED=y

cd /home/baohongg/afivo-streamer/programs/gas_gradient_3d/

# rerun
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="1.5 1 0 -4" -output%name=/scratch-shared/baohongg/3d_line/line1.5+x-4z_1.6e6_dr0.7_sw0.01 &
 
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="0.5 1 0 -2" -output%name=/scratch-shared/baohongg/3d_line/line0.5+x-2z_1.6e6_dr0.7_sw0.01 & 

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="0 1 0 -1" -output%name=/scratch-shared/baohongg/3d_line/linex-z_1.6e6_dr0.7_sw0.01 &

# new cases
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.03 -line_coeff="0.5 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.5-z_1.6e6_dr1.8_sw0.03 & 

wait

