#!/bin/env bash

# run_transform_pipeline.sh

# Define where the raw/blended files to which you want to register the atlas live 
raw_dir=$1
# Define the path to where elastix inverse transformations were run
elastix_atlas_to_auto_dir=$2
elastix_auto_to_cell_dir=$3
# Define where you want your aligned raw atlas files to live
output_dir=$4
# The raw-space atlas z planes will live in ${output_dir}/transformed_annotations/single_tifs/
annotation_volume_path=$5
dv=$6

# Set the atlas annotation volume that you want to warp to raw space


# STEP 1: copies TransformParameters.*.txt files, run transformix to get atlas into downsized space Single core 
OUT1=$(sbatch --parsable --export=ALL,raw_dir=${raw_dir},output_dir=${output_dir},elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir},annotation_volume_path=${annotation_volume_path} slurm_scripts/transform_step1.sh ) 
echo $OUT1

# STEP 2: Writes out z planes in raw space dimensions. Multi-core
OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,raw_dir=${raw_dir},elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir},output_dir=${output_dir},annotation_volume_path=${annotation_volume_path},dv=${dv} slurm_scripts/transform_step2.sh ) 
echo $OUT2
