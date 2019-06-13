#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 1440                 # time (minutes)
#SBATCH -o /scratch/zmd/logs/fullszdata_stack_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/fullszdata_stack_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 30000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

module load anacondapy/5.1.0
. activate lightsheet

python make_stacks_of_fullsizedata.py  #process zplns, check that 1000 > zplns/slurmfactor

# HOW TO USE:
# sbatch --array=0-20 sub_arrayjob.sh
#xvfb-run --auto-servernum --server-num=1