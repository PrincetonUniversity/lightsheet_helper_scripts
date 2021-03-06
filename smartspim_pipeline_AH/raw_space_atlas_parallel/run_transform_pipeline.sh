#!/bin/env bash

# run_transform_pipeline.sh

# Define where the raw/blended files to which you want to register the atlas live 
raw_dir=/jukebox/wang/Jess/lightsheet_output/202010_cfos/processed/an001/full_sizedatafld/072420_jv_ymazelearn_an1_1_3x_647_008na_1hfds_z10um_50msec_ch00
# Define the path to where elastix inverse transformations were run
elastix_atlas_to_auto_dir=/jukebox/wang/Jess/lightsheet_output/202010_cfos/processed/an001/ClearMapClusterOutput/elastix_auto_to_atlas
elastix_auto_to_cell_dir=/jukebox/wang/Jess/lightsheet_output/202010_cfos/processed/an001/ClearMapClusterOutput/elastix_cfos_to_auto
# Define where you want your aligned raw atlas files to live
output_dir=/jukebox/wang/ahoag/test_raw_atlas
# The raw-space atlas z planes will live in ${output_dir}/transformed_annotations/single_tifs/


# STEP 1: copies TransformParameters.*.txt files, run transformix to get atlas into downsized space Single core 
OUT1=$(sbatch --parsable --export=ALL,raw_dir=${raw_dir},output_dir=${output_dir},elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir} slurm_scripts/transform_step1.sh ) 
echo $OUT1

# STEP 2: Writes out z planes in raw space dimensions. Multi-core
OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,raw_dir=${raw_dir},elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir},output_dir=${output_dir}, slurm_scripts/transform_step2.sh ) 
echo $OUT2
