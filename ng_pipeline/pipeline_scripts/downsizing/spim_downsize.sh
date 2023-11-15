#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 6                      # number of cores
#SBATCH -t 150
#SBATCH -o logs/smartspim_downsize_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/smartspim_downsize_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 80000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

module load anacondapy/2020.11
. activate lightsheet

echo "Storage directory:" "${corrected_dir}"
echo "Destination directory:" "${downsized_dir}"

python spim_downsize.py "${corrected_dir}" "${downsized_dir}" "{dv}" "{atl}"
