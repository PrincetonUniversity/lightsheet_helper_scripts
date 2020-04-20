#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -n 1                      # number of cores
#SBATCH -t 60                 # time (minutes)
#SBATCH -o /scratch/zmd/logs/erode_%a.out        # STDOUT
#SBATCH -e /scratch/zmd/logs/erode_%a.err        # STDERR
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 20000

module load anacondapy/5.3.1
source activate lightsheet

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

echo "Array Index: $SLURM_ARRAY_TASK_ID"

python erode_all_structures.py ${SLURM_ARRAY_TASK_ID} 

# Usage notes:
# after = go once the specified job starts
# afterany = go if the specified job finishes, regardless of success
# afternotok = go if the specified job fails
# afterok = go if the specified job completes successfully
