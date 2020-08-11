#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 3                      # number of cores
#SBATCH -t 40
#SBATCH -o /scratch/zmd/logs/ts_thresh_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/ts_thresh_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 5000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

terastitcher --displthres --projin=/jukebox/LightSheetData/wang-mouse/seagravesk/20200810_13_10_58_f37080_mouse2_20171015_slow_focus/Ex_785_Em_3/xml_displproj.xml --projout=/jukebox/LightSheetData/wang-mouse/seagravesk/20200810_13_10_58_f37080_mouse2_20171015_slow_focus/Ex_785_Em_3/xml_displthres.xml --threshold=0.7
