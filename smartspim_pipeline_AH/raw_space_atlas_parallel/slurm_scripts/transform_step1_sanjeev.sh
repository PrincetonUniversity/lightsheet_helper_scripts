#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 60                # time (minutes)
#SBATCH -o logs/transform_step1_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/transform_step1_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 80000 #80 gbs

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

module load anacondapy/2020.11
module load elastix/4.8
. activate lightsheet

python transform_annotations_to_fullsize_cfos_sanjeev.py step1 ${raw_dir} \
	${elastix_atlas_to_auto_dir} ${elastix_auto_to_cell_dir} ${output_dir}
