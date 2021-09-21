#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 6                      # number of cores
#SBATCH -t 60
#SBATCH -o logs/downsize_%A_%a.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/downsize_%A_%a.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 # 100gbs

module load anacondapy/2020.11
. activate lightsheet

xvfb-run -d python downsizing/spim_downsize.py ${sample_dir} ${output_rootpath}
