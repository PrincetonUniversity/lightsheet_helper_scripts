#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 200                # time (minutes)
#SBATCH -o /scratch/kellyms/logs/convert_ann_%a_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e /scratch/kellyms/logs/convert_ann_%a_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 10000 #10 gbs

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

module load anacondapy/5.1.0
. activate lightsheet

echo "Array Index: $SLURM_ARRAY_TASK_ID"

python convert_memmap_npy_to_single_tifs.py ${SLURM_ARRAY_TASK_ID}
