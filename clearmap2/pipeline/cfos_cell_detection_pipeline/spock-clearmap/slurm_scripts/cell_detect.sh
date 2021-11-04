#!/bin/env bash
#
#BATCH -p all                # partition (queue)
#SBATCH -c 1                      # number of cores
#SBATCH -t 30                 # time (minutes)
#SBATCH -o logs/debug_cell_detect_%A_%a.out        # STDOUT #add _%a to see each array job
#SBATCH -e logs/debug_cell_detect_%A_%a.err        # STDERR #add _%a to see each array job
#SBATCH --contiguous #used to try and get cpu mem to be contigous
#SBATCH --mem 100000 #100gbs 

# tmpdir=/tmp/${UID}_clearmap2_block${SLURM_ARRAY_TASK_ID}
# echo "tmp dir is: $tmpdir"
authdir=./authfiles
export XDG_RUNTIME_DIR=$authdir

module load anacondapy/2020.11
conda activate ClearMap
# echo "Using display server: $displayport"
# Sleep for 2*array_job_id+1 seconds to avoid collisions between array jobs
tsleep=`echo "$SLURM_ARRAY_TASK_ID*2+1" | bc`
echo "Sleeping for $tsleep seconds"
sleep $tsleep
request_and_sample=`echo $sample_dir | cut -d "/" -f6,7`
cell_block_fname=${output_rootpath}/${request_and_sample}/${imaging_request}/rawdata/resolution_3.6x/cells_blocks/cells_block${SLURM_ARRAY_TASK_ID}.p
EXIT=0 # set in case file already exists. Then we skip while loop and just exit
while [[ ! -f "$cell_block_fname" ]] # use a while loop in case of BUS error when starting or finishing python script
do
	echo "In while loop, running cell detection on block ${SLURM_ARRAY_TASK_ID}"
	python spock-clearmap/clearmap_cell_detect.py ${sample_dir} ${imaging_request} ${blocks_per_job} ${output_rootpath} ${clearmap_params_file}
	EXIT=$?
done
echo "Exit status is $EXIT"
exit 0