#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 12                      # number of cores
#SBATCH -t 1000                # time (minutes)
#SBATCH -o /scratch/zmd/logs/kelly_reg_%a_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/zmd/logs/kelly_reg_%a_%j.err        # STDERR #add _%a to see each array job

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list
cat /proc/meminfo

module load anacondapy/5.3.1
. activate lightsheet

python kelly_downsize_and_reg.py ${SLURM_ARRAY_TASK_ID}