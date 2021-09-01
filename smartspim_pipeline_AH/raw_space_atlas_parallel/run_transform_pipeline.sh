#!/bin/env bash

# run_transform_pipeline.sh

# Define where the raw/blended files to which you want to register the atlas live 
raw_dir=/jukebox/LightSheetData/lightserv/oostland/MO_May2021_Tsc1_part1/MO_May2021_Tsc1_part1-521/imaging_request_1/output/processing_request_1/resolution_1.3x/full_sizedatafld/MO_521_3_AKA_20111303_flipped_ch01
# Define the path to where elastix inverse transformations were run
elastix_atlas_to_auto_dir=/jukebox/LightSheetData/lightserv/oostland/MO_May2021_Tsc1_part1/MO_May2021_Tsc1_part1-521/imaging_request_1/output/processing_request_1/resolution_1.3x/elastix_inverse_transform/injch_MO_521_3_AKA_20111303_flipped/MO_521_3_AKA_20111303_flipped_resized_ch00_resampledforelastix_atlas2reg
elastix_auto_to_cell_dir=/jukebox/LightSheetData/lightserv/oostland/MO_May2021_Tsc1_part1/MO_May2021_Tsc1_part1-521/imaging_request_1/output/processing_request_1/resolution_1.3x/elastix_inverse_transform/injch_MO_521_3_AKA_20111303_flipped/MO_521_3_AKA_20111303_flipped_resized_ch01_resampledforelastix_reg2sig
# Define where you want your aligned raw atlas files to live
output_dir=/jukebox/LightSheetData/lightserv/oostland/MO_May2021_Tsc1_part1/MO_May2021_Tsc1_part1-521/imaging_request_1/output/processing_request_1/resolution_1.3x/raw_atlas
# The raw-space atlas z planes will live in ${output_dir}/transformed_annotations/single_tifs/


# STEP 1: copies TransformParameters.*.txt files, run transformix to get atlas into downsized space Single core 
OUT1=$(sbatch --parsable --export=ALL,raw_dir=${raw_dir},output_dir=${output_dir},elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir} slurm_scripts/transform_step1.sh ) 
echo $OUT1

# STEP 2: Writes out z planes in raw space dimensions. Multi-core
OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,raw_dir=${raw_dir},elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir},output_dir=${output_dir}, slurm_scripts/transform_step2.sh ) 
echo $OUT2
