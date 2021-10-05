#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 6                  # number of cores
#SBATCH -t 200
#SBATCH -o logs/smartspim_inverse_reg_%A_%a.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/smartspim_inverse_reg_%A_%a.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 80000

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

module load anacondapy/2020.11
module load elastix/4.8
. activate lightsheet

xvfb-run -d python registration/spim_inverse_register.py ${sample_dir} ${imaging_request} ${output_rootpath}