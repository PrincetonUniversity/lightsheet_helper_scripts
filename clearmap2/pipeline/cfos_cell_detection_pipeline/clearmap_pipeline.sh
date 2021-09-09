#!/bin/env bash

#Stitching
echo "Experiment name / TeraStitcher folder hierarchy:" "$1"
echo "Storage directory:" "$2"

# # # import
OUT0=$(sbatch --parsable --export=ALL,input_dir=$1,output_dir=$2 stitching/slurm_scripts/ts_smartspim_import.sh)
echo $OUT0

# #displacement computation
OUT1=$(sbatch --parsable --dependency=afterok:${OUT0##* } --export=ALL,input_dir=$1,output_dir=$2 stitching/slurm_scripts/ts_smartspim_compute.sh)
echo $OUT1

# # #thresholding
OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,input_dir=$1,output_dir=$2 stitching/slurm_scripts/ts_smartspim_proj.sh)
echo $OUT2

#merge
OUT3=$(sbatch --parsable --dependency=afterok:${OUT2##* } --export=ALL,input_dir=$1,output_dir=$2 stitching/slurm_scripts/ts_smartspim_merge.sh)
echo $OUT3

#Pystripe
echo "Stitched channel path:" "$3"
echo "Flat path:" "$4"
echo "Corrected path:" "$5"

OUT4=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,stitched_dir=$3,flat_dir=$4,corrected_dir=$5 pystripe/pystripe.sh)
echo $OUT4

#ClearMap
echo "Sample path:" "$6"

#Preprocessing
OUT5=$(sbatch --parsable --dependency=afterok:${OUT4##* } --export=ALL,sample_dir=$6 spock-clearmap/slurm_scripts/preprocessing.sh)
echo $OUT5

#Cell Detection
OUT6=$(sbatch --parsable --dependency=afterok:${OUT5##* } --export=ALL,sample_dir=$6 spock-clearmap/slurm_scripts/main.sh)
echo $OUT6

#Downsizing
echo "Downsized path:" "$7"

OUT7=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,corrected_dir=$5,downsized_dir=$7 downsizing/spim_downsize.sh)

#Registration
echo "Registration channel:" "$8"
echo "Cell channel:" "$9"

#Inverse
#OUT8=$(sbatch --parsable --dependency=afterok:${OUT7##* } --export=ALL,src=$7,reg=$8,cell=$9 registration/slurm_scripts/spim_inverse_register.sh)
#echo $OUT8
#Normal
OUT8=$(sbatch --parsable --dependency=afterok:${OUT7##* } --export=ALL,src=$7,reg=$8,cell=$9 registration/slurm_scripts/spim_register.sh)
echo $OUT8

#ClearMap Postprocessing
OUT9=$(sbatch --parsable --dependency=afterok:${OUT6##* }:${OUT8##* } --export=ALL,sample_dir=$6 spock-clearmap/slurm_scripts/postprocessing.sh)
echo $OUT9


