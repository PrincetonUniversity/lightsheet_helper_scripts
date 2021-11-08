#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 2                      # number of cores
#SBATCH -t 1            # time (minutes)
#SBATCH -o logs/debug_clearmap_create_memmap_vol_%j.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/debug_clearmap_create_memmap_vol_%j.err        # STDERR #add _%a to see each array job

module load anacondapy/2020.11
conda activate ClearMap

request_and_sample=`echo $sample_dir | cut -d "/" -f6,7`
stitched_file_fname=${output_rootpath}/${request_and_sample}/${imaging_request}/rawdata/resolution_3.6x/stitched.npy
python spock-clearmap/clearmap_make_memmap_vol.py ${sample_dir} ${imaging_request} ${output_rootpath}
if [[ ! -f "$stitched_file_fname" ]]
then 
	echo "Stitched file not found after code ran. Error."
	exit 1
fi
exit 0