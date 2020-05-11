#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 12                      # number of cores
#SBATCH -t 1200                 # time (minutes)
#SBATCH -o /scratch/zmd/logs/whl_brn_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/whl_brn_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 #100gbs

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list
echo "Array Index: $SLURM_ARRAY_TASK_ID"

module load anacondapy/5.3.1
. activate lightsheet

python correct_clearmap_cell_detection.py ${SLURM_ARRAY_TASK_ID}
