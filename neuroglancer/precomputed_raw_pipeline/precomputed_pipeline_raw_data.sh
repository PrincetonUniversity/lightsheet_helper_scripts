#!/bin/env bash
#
# --- PURPOSE ---
# Pipeline to make precomputed (i.e. Neuroglancer-friendly) 
# volumes for raw light sheet images and
# detected cells in raw space

# author: Austin Hoag, mods by Z. D. 
# date: 04/08/2020

viz_dir=/jukebox/scratch/zmd/
animal_id=20170116_tp_bl6_lob6b_lpv_07

# First raw data
echo "Jobids for raw data precomputed steps 0, 1 and 2:"

#Step 0: Make info file and layer directory
OUT0=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_data_step0.sh) 
echo $OUT0

#Step 1: Upload raw data to vol (writes precomputed data to disk)

OUT1=$(sbatch --parsable --dependency=afterok:${OUT0##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_data_step1.sh) 
echo $OUT1

# Step 2: Make downsamples (higher mips) 
OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,viz_dir=${viz_dir},animal_id=${animal_id} precomputed_data_step2.sh) 
#echo $OUT2 