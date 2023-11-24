#!/bin/bash

#SBATCH --nodes=1	     # number of nodes on which to run
#SBATCH --ntasks=8	     # number of tasks to run
#SBATCH --cpus-per-task=16    # number of cpus required per task
#SBATCH --partition=rome	     # partition requested
#SBATCH --mem=128G             # memory requested
#SBATCH --time=5-00:00:00	     # wall-clock time limit


module load 2022
module load GCC/11.3.0

export OMP_NUM_THREADS=32
export GFORTRAN_UNBUFFERED_PRECONNECTED=y


cd /home/baohongg/humid_air/afivo-streamer/programs/humid_2d/

./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0000_fac1.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=10 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0010_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=10 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0010_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=10 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0010_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=10 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0010_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=11 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0011_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=11 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0011_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=11 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0011_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=11 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0011_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=12 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0012_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=12 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0012_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=12 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0012_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=12 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0012_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=13 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0013_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=13 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0013_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=13 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0013_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=13 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0013_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=14 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0014_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=14 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0014_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=14 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0014_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=14 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0014_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=15 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0015_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=15 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0015_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=15 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0015_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=15 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0015_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=16 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0016_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=16 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0016_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=16 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0016_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=16 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0016_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=17 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0017_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=17 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0017_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=17 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0017_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=17 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0017_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=18 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0018_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=18 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0018_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=18 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0018_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=18 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0018_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=19 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0019_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=19 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0019_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=19 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0019_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=19 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0019_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=20 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0020_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=20 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0020_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=20 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0020_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=20 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0020_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=21 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0021_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=21 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0021_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=21 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0021_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=21 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0021_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=22 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0022_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=22 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0022_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=22 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0022_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=22 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0022_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=23 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0023_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=23 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0023_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=23 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0023_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=23 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0023_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=24 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0024_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=24 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0024_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=24 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0024_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=24 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0024_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=25 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0025_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=25 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0025_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=25 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0025_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=25 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0025_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=26 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0026_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=26 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0026_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=26 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0026_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=26 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0026_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=27 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0027_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=27 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0027_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=27 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0027_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=27 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0027_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=28 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0028_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=28 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0028_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=28 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0028_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=28 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0028_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=29 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0029_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=29 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0029_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=29 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0029_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=29 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0029_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=30 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0030_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=30 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0030_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=30 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0030_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=30 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0030_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=31 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0031_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=31 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0031_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=31 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0031_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=31 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0031_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=32 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0032_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=32 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0032_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=32 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0032_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=32 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0032_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=33 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0033_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=33 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0033_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=33 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0033_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=33 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0033_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=34 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0034_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=34 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0034_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=34 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0034_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=34 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0034_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=35 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0035_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=35 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0035_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=35 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0035_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=35 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0035_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=36 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0036_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=36 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0036_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=36 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0036_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=36 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0036_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=37 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0037_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=37 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0037_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=37 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0037_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=37 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0037_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=38 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0038_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=38 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0038_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=38 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0038_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=38 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0038_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=39 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0039_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=39 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0039_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=39 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0039_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=39 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0039_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=40 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0040_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=40 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0040_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=40 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0040_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=40 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0040_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=41 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0041_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=41 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0041_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=41 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0041_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=41 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0041_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=42 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0042_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=42 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0042_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=42 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0042_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=42 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0042_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=43 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0043_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=43 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0043_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=43 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0043_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=43 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0043_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=44 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0044_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=44 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0044_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=44 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0044_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=44 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0044_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=45 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0045_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=45 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0045_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=45 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0045_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=45 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0045_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=46 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0046_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=46 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0046_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=46 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0046_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=46 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0046_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=47 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0047_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=47 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0047_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=47 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0047_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=47 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0047_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=48 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0048_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=48 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0048_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=48 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0048_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=48 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0048_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=49 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0049_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=49 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0049_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=49 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0049_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=49 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0049_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=50 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0050_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=50 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0050_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=50 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0050_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=50 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0050_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=51 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0051_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=51 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0051_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=51 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0051_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=51 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0051_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=52 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0052_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=52 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0052_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=52 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0052_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=52 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0052_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=53 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0053_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=53 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0053_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=53 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0053_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=53 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0053_fac0.75 &

wait

