#!/bin/env bash


src="$1"
reg="$2"
cell="$3"

# # normal registration
#OUT0=$(sbatch --parsable --export=ALL,src=${src},reg=${reg},cell=${cell} slurm_scripts/spim_register.sh)
#echo $OUT0

# inverse registration
 OUT1=$(sbatch --parsable --export=ALL,src=${src},reg=${reg},cell=${cell} slurm_scripts/spim_inverse_register.sh)
 echo $OUT1
# OUT1=$(sbatch --parsable --export=ALL,src=${src},reg=${reg} slurm_scripts/spim_inverse_register.sh)
# echo $OUT1
