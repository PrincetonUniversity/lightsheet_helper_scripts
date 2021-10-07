#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 25                 # time (minutes)
#SBATCH -o logs/clearmap_merge_blocks_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/clearmap_merge_blocks_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 #100gbs

module load anacondapy/2020.11
conda activate ClearMap

xvfb-run -d python spock-clearmap/clearmap_merge_blocks.py ${sample_dir} ${imaging_request} ${output_rootpath} ${clearmap_params_file}
