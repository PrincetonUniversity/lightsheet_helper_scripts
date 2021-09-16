#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 6                      # number of cores
#SBATCH -t 45                 # time (minutes)
#SBATCH -o logs/clearmap_postprocessing_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/clearmap_postprocessing_%j.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate ClearMap

xvfb-run -d python cz_clearmap_postprocessing.py ${sample_name}
