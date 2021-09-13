#!/bin/env bash

### How to run:
# clearmap_pipeline.sh ${corrected_dir} ${sample_name} ${downsized_dir} ${reg_channel} ${cell_channel} 

echo "Sample path:" "$2"

#Preprocessing
OUT5=$(sbatch --parsable --dependency=afterok:${OUT4##* } --export=ALL,sample_dir=$2 spock-clearmap/slurm_scripts/preprocessing.sh)
echo $OUT5

#Cell Detection
OUT6=$(sbatch --parsable --dependency=afterok:${OUT5##* } --export=ALL,sample_dir=$2 spock-clearmap/slurm_scripts/main.sh)
echo $OUT6

#Downsizing
echo "Downsized path:" "$3"

OUT7=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,corrected_dir=$1,downsized_dir=$3 downsizing/spim_downsize.sh)

#Registration
echo "Registration channel:" "$4"
echo "Cell channel:" "$5"

#Inverse
OUT8=$(sbatch --parsable --dependency=afterok:${OUT7##* } --export=ALL,src=$3,reg=$4,cell=$5 registration/slurm_scripts/spim_inverse_register.sh)
echo $OUT8
#Normal
# OUT8=$(sbatch --parsable --dependency=afterok:${OUT7##* } --export=ALL,src=$3,reg=$4,cell=$5 registration/slurm_scripts/spim_register.sh)
# echo $OUT8

#ClearMap Postprocessing
OUT9=$(sbatch --parsable --dependency=afterok:${OUT6##* }:${OUT8##* } --export=ALL,sample_dir=$2 spock-clearmap/slurm_scripts/postprocessing.sh)
echo $OUT9


