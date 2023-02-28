#!/bin/bash

#SBATCH --nodes=1	     # number of nodes on which to run
#SBATCH --ntasks=4	     # number of tasks to run
#SBATCH --cpus-per-task=32    # number of cpus required per task
#SBATCH --partition=thin	     # partition requested
#SBATCH --mem=120G             # memory requested
#SBATCH --time=5-00:00:00	     # wall-clock time limit


module load 2020
module load GCC/9.3.0

export OMP_NUM_THREADS=32
export GFORTRAN_UNBUFFERED_PRECONNECTED=y

cd afivo-streamer/programs/gas_gradient_3d/

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.02 -line_coeff="0.5 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.5-z_1.6e6_dr1.8_sw0.02 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.04 -line_coeff="0.5 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.5-z_1.6e6_dr1.8_sw0.04 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.06 -line_coeff="0.5 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.5-z_1.6e6_dr1.8_sw0.06 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.08 -line_coeff="0.5 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.5-z_1.6e6_dr1.8_sw0.08 &

wait

