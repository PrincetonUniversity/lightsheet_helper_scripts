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

terastitcher -1 --volin=/jukebox/LightSheetTransfer/kelly/20200630_14_03_08_m57207_dem_cfos_20190320/Ex_642_Em_2/ --ref1=1 --ref2=2 --ref3=3 --vxl1=1.866 --vxl2=1.866 --vxl3=4 --projout=xml_import
