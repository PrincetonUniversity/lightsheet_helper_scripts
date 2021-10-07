#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 25                 # time (minutes)
#SBATCH -o logs/diagnostic_merged_cells_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/diagnostic_merged_cells_%j.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate ClearMap

xvfb-run -d python spock-clearmap/diagnostic_merged_cells.py ${sample_dir} ${imaging_request} ${output_rootpath}