#!/bin/env bash

# echo variable names
atl="Allen"
image_resolution="3.6x"
dv="d"
step_size="2.0" 
channel_name="488"
channel_name2="647"
downsized_vol_name="downsized_for_atlas.tif"

viz_dir_488=$1
echo viz_dir_488
viz_dir_647=$2

corrected_dir_488=$3
downsized_dir_488=$4
corrected_dir_647=$5
downsized_dir_647=$6

src=$7

raw_dir=$8 #like raw files dir, different from src
output_dir=$9
elastix_atlas_to_auto_dir=${10}
elastix_auto_to_cell_dir=${11}

brain=${12}
raw_atlas_dir=${13}
viz_dir_ra=${14}



## Corrected Layers
echo "Make raw atlas layers"

# Make info file and layer directory
OUT0_c=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir_488},image_resolution=${image_resolution},\
channel_name=${channel_name} pipeline_scripts/corrected_pipeline/precomputed_step0.sh) 
echo $OUT0_c

OUT0_c2=$(sbatch --parsable --export=ALL,viz_dir=${viz_dir_647},image_resolution=${image_resolution},\
channel_name=${channel_name2} pipeline_scripts/corrected_pipeline/precomputed_step0.sh) 
echo $OUT0_c2

# Upload raw data to vol (writes precomputed data to disk)
OUT1_c=$(sbatch --parsable --dependency=afterok:${OUT0_c##* } --export=ALL,viz_dir=${viz_dir_488},\
image_resolution=${image_resolution},channel_name=${channel_name} \
pipeline_scripts/corrected_pipeline/precomputed_step1.sh) 
echo $OUT1_c

OUT1_c2=$(sbatch --parsable --dependency=afterok:${OUT0_c2##* } --export=ALL,viz_dir=${viz_dir_647},\
image_resolution=${image_resolution},channel_name=${channel_name2} \
pipeline_scripts/corrected_pipeline/precomputed_step1.sh) 
echo $OUT1_c2

# Transfer task -- rechunk the full resolution data 
OUT2_c=$(sbatch --parsable --dependency=afterok:${OUT1_c##* } --export=ALL,viz_dir=${viz_dir_488},\
image_resolution=${image_resolution},channel_name=${channel_name} \
pipeline_scripts/corrected_pipeline/precomputed_step2.sh) 
echo $OUT2_c

OUT2_c2=$(sbatch --parsable --dependency=afterok:${OUT1_c2##* } --export=ALL,viz_dir=${viz_dir_647},\
image_resolution=${image_resolution},channel_name=${channel_name2} \
pipeline_scripts/corrected_pipeline/precomputed_step2.sh) 
echo $OUT2_c2

# Downsampling 
OUT3_c=$(sbatch --parsable --dependency=afterok:${OUT2_c##* } --export=ALL,viz_dir=${viz_dir_488},\
image_resolution=${image_resolution},channel_name=${channel_name} \
pipeline_scripts/corrected_pipeline/precomputed_step3.sh) 
echo $OUT3_c

OUT3_c2=$(sbatch --parsable --dependency=afterok:${OUT2_c2##* } --export=ALL,viz_dir=${viz_dir_647},\
image_resolution=${image_resolution},channel_name=${channel_name2} \
pipeline_scripts/corrected_pipeline/precomputed_step3.sh) 
echo $OUT3_c2



## Downsiszing
# Downsize 488
OUT0=$(sbatch --parsable --export=ALL,corrected_dir_488=${corrected_dir_488},\
downsized_dir_488=${downsized_dir_488},dv=${dv},atl=${atl} \
pipeline_scripts/downsizing/spim_downsize.sh)
echo "Step 0: Downsizing 488 channel"
echo $OUT0

# Downsize 647
OUT1=$(sbatch --parsable --export=ALL,corrected_dir_647=${corrected_dir_647},\
downsized_dir_647=${downsized_dir_647},dv=${dv},atl=${atl} \
pipeline_scripts/downsizing/spim_downsize.sh)
echo "Step 1: Downsizing 647 channel"
echo $OUT0



## Inverse Registration
OUT2=$(sbatch --parsable --dependency=afterok:${OUT0##* }:${OUT1##* } \
--export=ALL,src=${src},reg=${downsized_dir_488},reg_vol=${downsized_vol_name},\
cell=${downsized_dir_647},cell_vol=${downsized_vol_name},atl=${atl} \
pipeline_scripts/registration/spim_inverse_register_semicustom.sh)
echo "Step 2: Inverse registration"
echo $OUT2



## Raw Atlas
# Raw Atlas Step 1
OUT3=$(sbatch --parsable --dependency=afterok:${OUT2##* } --export=ALL,raw_dir=${raw_dir},output_dir=${output_dir},\
elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir},\
annotation_volume_path=${atl} \
pipeline_scripts/raw_space_atlas_parallel/transform_step1.sh) 
echo "Step 3: Raw Space Atlas part 1"
echo $OUT3

# Raw Atlas Step 2
OUT4=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,raw_dir=${raw_dir},\
elastix_atlas_to_auto_dir=${elastix_atlas_to_auto_dir},elastix_auto_to_cell_dir=${elastix_auto_to_cell_dir},\
output_dir=${output_dir},annotation_volume_path=${atl},dv=${dv} \
pipeline_scripts/raw_space_atlas_parallel/transform_step2.sh) 
echo "Step 4: Raw Space Atlas part 2"
echo $OUT4



## Make Raw Atlas Layers
echo "Make raw atlas layers"
# Step 0: Make info file and layer directory
OUT5=$(sbatch --parsable --dependency=afterok:${OUT4##* } --export=ALL,brain=${brain},raw_atlas_dir=${raw_atlas_dir},\
viz_dir=${viz_dir_ra},step_size=${step_size},image_resolution=${image_resolution} \
pipeline_scripts/raw_atlas_pipeline_sanjeev/precomputed_atlas_step0.sh) 
echo $OUT5

# # Step 1: Upload raw data to vol (writes precomputed data to disk)
OUT6=$(sbatch --parsable --dependency=afterok:${OUT5##* } --export=ALL,brain=${brain},\
raw_atlas_dir=${raw_atlas_dir},viz_dir=${viz_dir_ra},step_size=${step_size},\
image_resolution=${image_resolution} pipeline_scripts/raw_atlas_pipeline_sanjeev/precomputed_atlas_step1.sh) 
echo $OUT6

#Step 2: Transfer tasks 
OUT7=$(sbatch --parsable --dependency=afterok:${OUT6##* } --export=ALL,brain=${brain},\
raw_atlas_dir=${raw_atlas_dir},viz_dir=${viz_dir_ra},step_size=${step_size},\
image_resolution=${image_resolution} pipeline_scripts/raw_atlas_pipeline_sanjeev/precomputed_atlas_step2.sh) 
echo $OUT7

#Step 3: Make downsamples (higher mips) 
OUT8=$(sbatch --parsable --dependency=afterok:${OUT7##* } --export=ALL,brain=${brain},\
raw_atlas_dir=${raw_atlas_dir},viz_dir=${viz_dir_ra},step_size=${step_size},\
image_resolution=${image_resolution} pipeline_scripts/raw_atlas_pipeline_sanjeev/precomputed_atlas_step3.sh) 
echo $OUT8


