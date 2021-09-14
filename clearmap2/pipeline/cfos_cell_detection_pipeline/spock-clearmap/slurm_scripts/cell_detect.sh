#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 1                 # time (minutes)
#SBATCH -o logs/clearmap_cell_detect_%A_%a.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/clearmap_cell_detect_%A_%a.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate ClearMap

xvfb-run -d python cz_clearmap_cell_detect.py ${sample_name}
