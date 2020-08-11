#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 3                      # number of cores
#SBATCH -t 800
#SBATCH -o /scratch/zmd/logs/ts_merge_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/ts_merge_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 85000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

terastitcher --merge --projin=/jukebox/LightSheetData/wang-mouse/seagravesk/20200810_13_10_58_f37080_mouse2_20171015_slow_focus/Ex_785_Em_3/xml_placetiles.xml --volout=/jukebox/LightSheetData/wang-mouse/seagravesk/20200810_13_10_58_f37080_mouse2_20171015_slow_focus/Ex_785_Em_3/stitched --imout_depth=16 --resolutions=0
