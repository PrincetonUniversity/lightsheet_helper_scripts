#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 3                      # number of cores
#SBATCH -t 10
#SBATCH -o /scratch/zmd/logs/ts_import_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/ts_import_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 5000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

terastitcher -1 --volin=/jukebox/LightSheetData/wang-mouse/seagravesk/20200810_13_10_58_f37080_mouse2_20171015_slow_focus/Ex_785_Em_3 --ref1=x --ref2=y --ref3=z --vxl1=1.81 --vxl2=1.81 --vxl3=2 --projout=xml_import
