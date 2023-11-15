#!/bin/env bash
#
# --- PURPOSE ---
# Pipeline to make precomputed (i.e. Neuroglancer-friendly) 
# volumes for raw atlas 

# author: Austin Hoag
# date: 04/08/2020

brain=$1
raw_atlas_dir=$2
viz_dir=$3
step_size=$4
image_resolution=$5
# dataset=201810_adultacutePC_ymaze_cfos
# animal_id=dadult_pc_crus1_5
# dataset=202002_cfos
# animal_id=an5_cno
# dataset=201904_ymaze_cfos
# animal_id=an27 


# echo "Jobids for raw atlas precomputed steps 0, 1 and 2:"

# Step 0: Make info file and layer directory
OUT0=$(sbatch --parsable --export=ALL,brain=${brain},raw_atlas_dir=${raw_atlas_dir},viz_dir=${viz_dir},step_size=${step_size},image_resolution=${image_resolution} precomputed_atlas_step0.sh) 
echo $OUT0

# # Step 1: Upload raw data to vol (writes precomputed data to disk)
OUT1=$(sbatch --parsable --dependency=afterok:${OUT0##* } --export=ALL,brain=${brain},raw_atlas_dir=${raw_atlas_dir},viz_dir=${viz_dir},step_size=${step_size},image_resolution=${image_resolution} precomputed_atlas_step1.sh) 
echo $OUT1

#Step 2: Transfer tasks 
OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,brain=${brain},raw_atlas_dir=${raw_atlas_dir},viz_dir=${viz_dir},step_size=${step_size},image_resolution=${image_resolution} precomputed_atlas_step2.sh) 
echo $OUT2

#Step 3: Make downsamples (higher mips) 
OUT3=$(sbatch --parsable --dependency=afterok:${OUT2##* } --export=ALL,brain=${brain},raw_atlas_dir=${raw_atlas_dir},viz_dir=${viz_dir},step_size=${step_size},image_resolution=${image_resolution} precomputed_atlas_step3.sh) 
echo $OUT3

# Usage notes:
# after = go once the specified job starts
# afterany = go if the specified job finishes, regardless of success
# afternotok = go if the specified job fails
# afterok = go if the specified job completes successfully
