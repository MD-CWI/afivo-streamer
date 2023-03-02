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

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.9 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/scratch-shared/baohongg/3d_sphere/in_1.6e6_dr0.9_sw0.01_sc0.7-0.5_sr0.2 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/scratch-shared/baohongg/3d_sphere/in_1.6e6_dr0.8_sw0.01_sc0.7-0.5_sr0.2 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/scratch-shared/baohongg/3d_sphere/in_1.6e6_dr0.7_sw0.01_sc0.7-0.5_sr0.2 &

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/scratch-shared/baohongg/3d_sphere/in_1.6e6_dr0.6_sw0.01_sc0.7-0.5_sr0.2 &

wait

