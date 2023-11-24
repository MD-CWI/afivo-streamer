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

./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=54 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0054_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=54 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0054_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=54 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0054_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=54 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0054_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=55 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0055_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=55 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0055_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=55 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0055_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=55 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0055_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=56 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0056_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=56 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0056_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=56 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0056_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=56 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0056_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=57 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0057_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=57 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0057_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=57 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0057_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=57 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0057_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=58 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0058_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=58 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0058_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=58 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0058_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=58 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0058_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=59 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0059_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=59 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0059_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=59 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0059_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=59 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0059_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=60 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0060_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=60 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0060_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=60 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0060_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=60 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0060_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=61 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0061_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=61 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0061_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=61 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0061_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=61 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0061_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=62 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0062_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=62 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0062_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=62 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0062_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=62 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0062_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=63 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0063_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=63 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0063_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=63 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0063_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=63 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0063_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=64 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0064_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=64 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0064_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=64 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0064_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=64 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0064_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=65 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0065_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=65 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0065_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=65 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0065_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=65 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0065_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=66 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0066_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=66 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0066_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=66 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0066_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=66 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0066_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=67 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0067_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=67 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0067_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=67 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0067_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=67 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0067_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=68 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0068_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=68 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0068_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=68 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0068_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=68 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0068_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=69 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0069_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=69 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0069_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=69 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0069_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=69 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0069_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=70 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0070_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=70 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0070_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=70 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0070_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=70 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0070_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=71 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0071_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=71 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0071_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=71 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0071_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=71 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0071_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=72 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0072_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=72 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0072_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=72 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0072_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=72 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0072_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=73 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0073_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=73 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0073_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=73 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0073_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=73 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0073_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=74 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0074_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=74 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0074_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=74 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0074_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=74 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0074_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=75 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0075_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=75 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0075_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=75 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0075_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=75 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0075_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=76 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0076_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=76 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0076_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=76 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0076_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=76 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0076_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=77 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0077_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=77 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0077_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=77 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0077_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=77 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0077_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=78 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0078_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=78 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0078_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=78 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0078_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=78 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0078_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=79 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0079_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=79 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0079_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=79 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0079_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=79 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0079_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=80 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0080_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=80 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0080_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=80 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0080_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=80 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0080_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=81 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0081_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=81 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0081_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=81 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0081_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=81 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0081_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=82 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0082_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=82 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0082_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=82 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0082_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=82 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0082_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=83 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0083_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=83 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0083_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=83 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0083_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=83 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0083_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=84 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0084_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=84 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0084_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=84 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0084_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=84 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0084_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=85 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0085_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=85 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0085_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=85 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0085_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=85 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0085_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=86 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0086_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=86 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0086_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=86 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0086_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=86 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0086_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=87 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0087_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=87 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0087_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=87 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0087_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=87 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0087_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=88 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0088_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=88 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0088_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=88 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0088_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=88 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0088_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=89 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0089_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=89 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0089_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=89 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0089_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=89 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0089_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=90 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0090_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=90 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0090_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=90 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0090_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=90 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0090_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=91 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0091_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=91 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0091_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=91 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0091_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=91 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0091_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=92 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0092_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=92 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0092_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=92 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0092_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=92 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0092_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=93 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0093_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=93 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0093_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=93 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0093_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=93 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0093_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=94 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0094_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=94 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0094_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=94 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0094_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=94 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0094_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=95 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0095_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=95 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0095_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=95 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0095_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=95 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0095_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=96 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0096_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=96 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0096_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=96 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0096_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=96 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0096_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=97 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0097_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=97 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0097_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=97 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0097_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=97 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0097_fac0.75 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=98 -input_data%modified_rate_factors=0.0 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0098_fac0.0 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=98 -input_data%modified_rate_factors=0.25 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0098_fac0.25 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=98 -input_data%modified_rate_factors=0.5 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0098_fac0.5 &
./streamer d200mm_sensitivity_malagon_0.03H2O.cfg -input_data%modified_reaction_ix=98 -input_data%modified_rate_factors=0.75 -output%name=/scratch-shared/baohongg/malagon_chemistry_ix0098_fac0.75 &

wait

