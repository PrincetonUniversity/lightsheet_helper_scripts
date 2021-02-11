#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 10
#SBATCH -o logs/ts_proj_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/ts_proj_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous

# echo "In the directory: `pwd` "
# echo "As the user: `whoami` "
# echo "on host: `hostname` "

# cat /proc/$$/status | grep Cpus_allowed_list

module load anacondapy/5.3.1
. activate lightsheet

xvfb-run -d python spim_stitch.py step2 ${input_dir} ${output_dir} 
