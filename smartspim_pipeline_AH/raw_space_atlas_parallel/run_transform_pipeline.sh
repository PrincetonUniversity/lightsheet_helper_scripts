#!/bin/env bash

# run_transform_pipeline.sh -- a bash script for create a raw-space atlas

# Define where the raw/blended files to which you want to register the atlas live 
raw_dir=/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_14_20_11_f37073_mouse1_20171010/Ex_785_Em_3/corrected

# Define the path to where elastix inverse transformations were run
elastix_atlas_to_auto_dir=/jukebox/wang/seagravesk/lightsheet/201710_cfos_left_side_only_registration/f37073_mouse1/elastix_inverse_transform/cellch_f37073_mouse1_20171010_790_015na_1hfsds_z5um_1000msec/f37073_mouse1_20171010_488_015na_1hfsds_z5um_150msec_resized_ch00_resampledforelastix_atlas2reg
elastix_auto_to_cell_dir=/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_14_20_11_f37073_mouse1_20171010/elastix_inverse_transform/

# Define where you want your aligned raw atlas files to live
# The raw-space atlas z planes will live in ${output_dir}/transformed_annotations/single_tifs/
output_dir=/jukebox/LightSheetData/wang-mouse/seagravesk/20200901_14_20_11_f37073_mouse1_20171010/Ex_785_Em_3/raw_atlas

# Set the atlas annotation volume that you want to warp to raw space
#annotation_volume_path=/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing.tif
annotation_volume_path=/jukebox/LightSheetData/atlas/annotation_sagittal_atlas_20um_16bit_hierarch_labels.tif

# STEP 1: copies TransformParameters.*.txt files, run transformix to get atlas into downsized space Single core 
OUT1=$(sbatch --parsable \
	--export=ALL,raw_dir=${raw_dir},output_dir=${output_dir},elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir},annotation_volume_path=${annotation_volume_path} \
	slurm_scripts/transform_step1.sh ) 
echo $OUT1

# STEP 2: Writes out z planes in raw space dimensions. Multi-core
OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } \
	--export=ALL,raw_dir=${raw_dir},output_dir=${output_dir},elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir},annotation_volume_path=${annotation_volume_path} \
	 slurm_scripts/transform_step2.sh ) 
echo $OUT2
