#!/bin/env bash 
# 
#SBATCH -p all # partition (queue) 
#SBATCH -c 1 # number of cores 
#SBATCH -t 200 # time (minutes) 
#SBATCH -o /jukebox/scratch/zmd/logs/cellcnt_%j.out # STDOUT #add _%a to see each array job 
#SBATCH -e /jukebox/scratch/zmd/logs/cellcnt_%j.err # STDERR #add _%a to see each array job 
#SBATCH --contiguous #used to try and get cpu mem to be contigous 
#SBATCH --mem 50000 #50 gbs 

module load anacondapy/5.1.0 
. activate idisco

python 20190415_jv_analyze_cell_counts.py
