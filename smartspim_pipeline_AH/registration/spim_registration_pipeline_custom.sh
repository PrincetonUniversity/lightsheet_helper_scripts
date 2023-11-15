#!/bin/env bash


src="$1"
atl="$2"
out="$3"

# # normal registration
OUT0=$(sbatch --parsable --export=ALL,src=${src},atl=${atl},out=${out} slurm_scripts/spim_register_custom.sh)
echo $OUT0

