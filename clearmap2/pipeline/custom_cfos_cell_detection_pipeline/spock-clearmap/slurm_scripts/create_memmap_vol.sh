#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 12                      # number of cores
#SBATCH -t 150                 # time (minutes)
#SBATCH -o logs/clearmap_create_memmap_vol_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/clearmap_create_memmap_vol_%j.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate ClearMap

sample=$(basename ${sample_dir})
stitched_file_fname=${output_rootpath}/${request_name}/${sample}/stitched.npy
python spock-clearmap/clearmap_make_memmap_vol.py ${sample_dir} ${request_name} ${output_rootpath}
echo $stitched_file_fname
if [[ ! -f "$stitched_file_fname" ]]
then 
	echo "Stitched file not found after code ran. Error."
	exit 1
fi
exit 0