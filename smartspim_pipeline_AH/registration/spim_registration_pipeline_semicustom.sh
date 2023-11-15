#!/bin/env bash


src="$1"
reg="$2"
reg_vol="$3"
cell="$4"
cell_vol="$5"
atl="$6"
dir="$7"
#echo $dir

# # normal registration
if [[ "$dir" == "forward" ]]
then
	OUT0=$(sbatch --parsable --export=ALL,src=${src},reg=${reg},reg_vol=${reg_vol},cell=${cell},cell_vol=${cell_vol},atl=${atl} slurm_scripts/spim_register_semicustom.sh)
	echo $OUT0
fi

# inverse registration
if [[ "$dir" == "inverse" ]]
then
	OUT1=$(sbatch --parsable --export=ALL,src=${src},reg=${reg},reg_vol=${reg_vol},cell=${cell},cell_vol=${cell_vol},atl=${atl} slurm_scripts/spim_inverse_register_semicustom.sh)
	echo $OUT1
fi

