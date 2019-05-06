#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 6                      # number of cores
#SBATCH -t 1440                 # time (minutes)
#SBATCH -o /scratch/zmd/logs/step1_terastitcher_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/step1_terastitcher_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 80000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

module load anacondapy/2.7
. activate lightsheet

xvfb-run python terastitch.py  #process zplns, check that 1000 > zplns/slurmfactor

# HOW TO USE:
# sbatch --array=0-20 sub_arrayjob.sh
#xvfb-run --auto-servernum --server-num=1
