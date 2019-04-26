#!/bin/env bash 
# 
#SBATCH -p all # partition (queue) 
#SBATCH -c 1 # number of cores 
#SBATCH -t 200 # time (minutes) 
#SBATCH -o /jukebox/scratch/zmd/logs/qc_%j.out # STDOUT #add _%a to see each array job 
#SBATCH -e /jukebox/scratch/zmd/logs/qc_%j.err # STDERR #add _%a to see each array job 
#SBATCH --contiguous #used to try and get cpu mem to be contigous 
#SBATCH --mem 30000 #30 gbs 

module load anacondapy/5.3.1 
. activate lspy3 

python 20190422_jv_qc.py
