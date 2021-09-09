#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 8                 # number of cores
#SBATCH -t 150                 # number of minutes 
#SBATCH -o logs/spim_pystripe_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/spim_pystripe_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 25000                      #RAM (MBs)- 25GBS

cat /proc/$$/status | grep Cpus_allowed_list

#required
module load anacondapy/2020.11
conda activate lightsheet

pystripe -i ${stitched_dir} -f ${flat_dir} -o ${corrected_dir}
