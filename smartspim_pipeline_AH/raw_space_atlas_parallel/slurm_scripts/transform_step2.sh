#!/bin/env bash
#
#SBATCH -p all                # partition (queue)
#SBATCH -c 4                      # number of cores
#SBATCH -t 150               # time (minutes)
#SBATCH -o logs/transform_step2_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/transform_step2_%j.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 #100 gbs

echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "

cat /proc/$$/status | grep Cpus_allowed_list

module load anacondapy/2020.11
module load elastix/4.8
. activate lightsheet

python transform_annotations_to_fullsize_cfos.py step2 ${raw_dir} \
	${elastix_atlas_to_auto_dir} ${elastix_auto_to_cell_dir} ${output_dir} \
	${annotation_volume_path} ${dv}
