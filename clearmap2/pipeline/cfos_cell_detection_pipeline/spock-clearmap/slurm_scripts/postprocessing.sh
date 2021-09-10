#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 12                      # number of cores
#SBATCH -t 1200                 # time (minutes)
#SBATCH -o /scratch/zmd/logs/whl_brn_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/whl_brn_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 #100gbs

module load anacondapy/2020.11
conda activate clearmap

xvfb-run -d python cz_clearmap_postprocessing.py ${sample_dir}
