#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -n 1                      # number of cores
#SBATCH -t 60                 # time (minutes)
#SBATCH -o /scratch/zmd/logs/erode_sum_v2.out        # STDOUT
#SBATCH -e /scratch/zmd/logs/erode_sum_v2.err        # STDERR
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 50000

module load anacondapy/5.3.1
source activate lightsheet

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

python sum_eroded_annotation_files.py

# Usage notes:
# after = go once the specified job starts
# afterany = go if the specified job finishes, regardless of success
# afternotok = go if the specified job fails
# afterok = go if the specified job completes successfully
