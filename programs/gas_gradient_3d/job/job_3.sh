#!/bin/bash

#SBATCH --nodes=1	     # number of nodes on which to run
#SBATCH --ntasks=4	     # number of tasks to run
#SBATCH --cpus-per-task=32    # number of cpus required per task
#SBATCH --partition=thin	     # partition requested
#SBATCH --mem=120G             # memory requested
#SBATCH --time=5-00:00:00	     # wall-clock time limit


module load 2022
module load GCC/11.3.0

export OMP_NUM_THREADS=32
export GFORTRAN_UNBUFFERED_PRECONNECTED=y

cd /home/baohongg/afivo-streamer/programs/gas_gradient_3d/

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=6 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.35-z_1.6e6_dr6_sw0.35 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=8 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.35-z_1.6e6_dr8_sw0.35 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=10 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.35-z_1.6e6_dr10_sw0.35 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=20 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/scratch-shared/baohongg/3d_line/line0.35-z_1.6e6_dr20_sw0.35 &

wait

