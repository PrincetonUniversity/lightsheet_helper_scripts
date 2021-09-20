#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 1                 # number of cores
#SBATCH -t 5                 # time (minutes)
#SBATCH -o logs/diagnostic_downsize_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/diagnostic_downsize_%j.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate ClearMap

xvfb-run -d python downsizing/diagnostic_downsize.py ${sample_dir} ${output_rootpath}
