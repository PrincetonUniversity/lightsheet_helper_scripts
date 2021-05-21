#!/bin/env bash

# Run this file like: 
# spim_stitching_pipeline.sh /path/to/terastitcher_folder_hierarchy /path/to/stitched

echo "Experiment name / TeraStitcher folder hierarchy:" "$1"
echo "Storage directory:" "$2"

# # # import
# OUT0=$(sbatch --parsable --export=ALL,input_dir=$1,output_dir=$2 slurm_scripts/ts_smartspim_import.sh)
# echo $OUT0

# # #displacement computation
# OUT1=$(sbatch --parsable --dependency=afterok:${OUT0##* } --export=ALL,input_dir=$1,output_dir=$2 slurm_scripts/ts_smartspim_compute.sh)
# echo $OUT1

# # # #thresholding
# OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,input_dir=$1,output_dir=$2 slurm_scripts/ts_smartspim_proj.sh)
# echo $OUT2

# #merge
# OUT3=$(sbatch --parsable --dependency=afterok:${OUT2##* } --export=ALL,input_dir=$1,output_dir=$2 slurm_scripts/ts_smartspim_merge.sh)
# echo $OUT3
OUT3=$(sbatch --parsable --export=ALL,input_dir=$1,output_dir=$2 slurm_scripts/ts_smartspim_merge.sh)
echo $OUT3



