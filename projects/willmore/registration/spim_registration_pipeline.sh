#!/bin/env bash

src='/jukebox/LightSheetData/lightserv/soline/2020.10.05_2357_819/2020.10.05_2357_819-001/imaging_request_1/rawdata/resolution_3.6x/20210209_16_17_16_M2357/'
reg='Ex_488_Em_0_corrected'
cell='Ex_561_Em_1_corrected'

# # normal registration
# OUT0=$(sbatch --parsable --export=ALL,src=${src},reg=${reg},cell=${cell} slurm_scripts/spim_register.sh)
# echo $OUT0

# inverse registration
# OUT1=$(sbatch --parsable --export=ALL,src=${src},reg=${reg},cell=${cell} slurm_scripts/spim_inverse_register.sh)
# echo $OUT1
OUT1=$(sbatch --parsable --export=ALL,src=${src},reg=${reg} slurm_scripts/spim_inverse_register.sh)
echo $OUT1