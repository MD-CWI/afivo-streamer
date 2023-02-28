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

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.6 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/scratch-shared/baohongg/3d_sphere/in_1.6e6_dr0.7_sw0.01_sc0.6-0.5_sr0.2 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/scratch-shared/baohongg/3d_sphere/in_1.6e6_dr0.7_sw0.01_sc0.5-0.5_sr0.2 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.2 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/scratch-shared/baohongg/3d_sphere/out_1.6e6_dr1.2_sw0.01_sc0.5-1.0_sr0.5 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/scratch-shared/baohongg/3d_sphere/out_1.6e6_dr1.4_sw0.01_sc0.5-1.0_sr0.5 &

wait

