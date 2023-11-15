#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH --nodes=1             # node count
#SBATCH -n 12                 # number of cores
#SBATCH -t 300                 # time (minutes)
#SBATCH -o logs/precomputed_step3_%j.out        # STDOUT
#SBATCH -e logs/precomputed_step3_%j.err        # STDERR


# start=$(date +%s.%N)
# echo "In the directory: `pwd` "
# echo "As the user: `whoami` "
# echo "on host: `hostname` "

# cat /proc/$$/status | grep Cpus_allowed_list
# cat /proc/meminfo

# echo "Array Allocation Number: $SLURM_ARRAY_JOB_ID"
# echo "Array Index: $SLURM_ARRAY_TASK_ID"

module load anacondapy/2020.11
. activate precomputed
xvfb-run -d python make_precomputed_corrected.py step3 ${viz_dir} ${image_resolution} ${channel_name} 

# finish=$(date +%s.%N)
# echo "$finish $start" | awk '{print "took " $1-$2 " seconds"}'

# HOW TO USE:
# sbatch --array=0-33 precomputed_step1.sh

# Usage notes:
# after = go once the specified job starts
# afterany = go if the specified job finishes, regardless of success
# afternotok = go if the specified job fails
# afterok = go if the specified job completes successfully
