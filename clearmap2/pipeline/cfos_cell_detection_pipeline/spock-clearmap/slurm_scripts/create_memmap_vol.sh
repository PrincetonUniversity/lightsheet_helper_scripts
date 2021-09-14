#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 8                      # number of cores
#SBATCH -t 100                 # time (minutes)
#SBATCH -o logs/clearmap_preprocessing_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/clearmap_preprocessing_%j.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate ClearMap

xvfb-run -d python spock-clearmap/cz_clearmap_make_memmap.py ${sample_dir}
