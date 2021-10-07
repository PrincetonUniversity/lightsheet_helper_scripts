#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 8                      # number of cores
#SBATCH -t 60                 # time (minutes)
#SBATCH -o logs/clearmap_create_memmap_vol_test_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/clearmap_create_memmap_vol_test_%j.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate ClearMap

xvfb-run -d python spock-clearmap/clearmap_make_memmap_vol.py ${sample_dir} ${imaging_request} ${output_rootpath}