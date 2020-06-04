#!/bin/env bash
#
# --- PURPOSE ---
# Pipeline to make precomputed (i.e. Neuroglancer-friendly) 
# volumes for raw light sheet images and
# detected cells in raw space

# author: Austin Hoag
# date: 04/08/2020

viz_dir=/jukebox/LightSheetData/lightserv_testing/neuroglancer/201904_ymaze_cfos/
animal_id=31

# First raw data
echo "Jobids for raw data precomputed steps 0, 1 and 2:"

# #Step 0: Make info file and layer directory
# OUT0=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_data_step0.sh) 
# echo $OUT0

# #Step 1: Upload raw data to vol (writes precomputed data to disk)

# OUT1=$(sbatch --parsable --dependency=afterok:${OUT0##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_data_step1.sh) 
# echo $OUT1

# # Step 2: Make downsamples (higher mips) 
# OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_data_step2.sh) 
# echo $OUT2 

# ## Now raw cells
# echo "Jobids for raw cells precomputed steps 0, 1 and 2:"
# #Step 0: Make info file and layer directory

# OUT3=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_cells_step0.sh) 
# echo $OUT3

# #Step 1: Upload raw cells to vol (writes precomputed data to disk)
# OUT4=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_cells_step1.sh) 
# echo $OUT4

# Step 2: Make downsamples (higher mips) 
# OUT5=$(sbatch --parsable --dependency=afterok:${OUT4##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_cells_step2.sh) 
# echo $OUT5
OUT5=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_cells_step2.sh) 
echo $OUT5

# ## Now raw atlas
# echo "Jobids for raw atlas precomputed steps 0, 1 and 2:"

# #Step 0: Make info file and layer directory
# OUT6=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_atlas_step0.sh) 
# echo $OUT6

# #Step 1: Upload raw data to vol (writes precomputed data to disk)
# OUT7=$(sbatch --parsable --dependency=afterok:${OUT6##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_atlas_step1.sh) 
# echo $OUT7

# #Step 2: Make downsamples (higher mips) 
# OUT8=$(sbatch --parsable --dependency=afterok:${OUT7##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_atlas_step2.sh) 
# echo $OUT8

# Usage notes:
# after = go once the specified job starts
# afterany = go if the specified job finishes, regardless of success
# afternotok = go if the specified job fails
# afterok = go if the specified job completes successfully
