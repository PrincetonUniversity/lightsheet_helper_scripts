#!/bin/env bash

### How to run:
# clearmap_pipeline.sh ${corrected_dir} ${sample_name} ${downsized_dir} ${reg_channel} ${cell_channel} 

echo "Corrected dir: $1"
echo "Sample name: $2"
echo "Downsized dir: $3"
echo "Registration channel: $4"
echo "Cell channel: $5"

#Preprocessing
OUT1=$(sbatch --parsable --dependency=afterok:${OUT4##* } --export=ALL,sample_name=$2 spock-clearmap/slurm_scripts/preprocessing.sh)
echo $OUT1

#Cell Detection
OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,sample_name=$2 spock-clearmap/slurm_scripts/main.sh)
echo $OUT2

#Downsizing
OUT3=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,corrected_dir=$1,downsized_dir=$3 downsizing/spim_downsize.sh)
echo $OUT3

#Inverse registration
OUT4=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,src=$3,reg=$4,cell=$5 registration/slurm_scripts/spim_inverse_register.sh)
echo $OUT4

#ClearMap Postprocessing
OUT5=$(sbatch --parsable --dependency=afterok:${OUT2##* }:${OUT4##* } --export=ALL,sample_name=$2 spock-clearmap/slurm_scripts/postprocessing.sh)
echo $OUT5

