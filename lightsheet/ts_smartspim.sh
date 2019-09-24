#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 3                      # number of cores
#SBATCH -t 400
#SBATCH -o logs/ts_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/ts_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mail-type=begin        # send mail when process begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=zmd@princeton.edu
#SBATCH --mem 35000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

echo "Array Allocation Number: $SLURM_ARRAY_JOB_ID"
echo "Array Index: $SLURM_ARRAY_TASK_ID"

terastitcher --merge --projin='/jukebox/LightSheetTransfer/microscope_tests/20190905_19_17_27_Princeton-4x_tiffs/Ex_642/xml_placetiles.xml' --volout='/jukebox/LightSheetTransfer/microscope_tests/20190905_19_17_27_Princeton-4x_tiffs/Ex_642/ts_out' --imout_depth=16 --resolutions=012345

# HOW TO USE:
# sbatch --array=0-20 sub_arrayjob.sh 
#xvfb-run --auto-servernum --server-num=1 
