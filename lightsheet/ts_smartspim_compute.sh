#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 3                      # number of cores
#SBATCH -t 800
#SBATCH -o /scratch/zmd/logs/ts_compute_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/ts_compute_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 85000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

terastitcher --displcompute --projin=/jukebox/LightSheetTransfer/kelly/20200630_14_03_08_m57207_dem_cfos_20190320/Ex_642_Em_2/xml_import.xml --algorithm=MIPNCC --sV=25 --sH=25 --sD=5

