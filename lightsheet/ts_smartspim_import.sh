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

terastitcher -1 --volin=/jukebox/LightSheetData/brodyatlas/raw_data/200728_k315_1_1x_488_016na_1hfds_z10um_50msec_20povlp_16-34-15/k315_1_1x_488_016na_1hfds_z10um_50msec_20povlp_ch00_left_lightsheet --ref1=x --ref2=y --ref3=z --vxl1=5.91 --vxl2=5.91 --vxl3=10 --projout=xml_import
