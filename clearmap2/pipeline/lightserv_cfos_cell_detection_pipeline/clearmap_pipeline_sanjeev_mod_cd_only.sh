#!/bin/env bash

### clearmap_pipeline.sh -- detect cells, register them to an atlas
# and create a dataframe of counts in brain atlas regions.
# This script is designed to process all samples in a single imaging request of a given request

### How to run:
# 
# clearmap_pipeline.sh $request_name $imaging_request $sample_name
# 
# The sample_name argument is optional. Include it if you want to run the pipeline ONLY for that one sample
# If you omit sample_name, then the pipeline will run run all samples in a request at a given imaging request
# For example if you want to run for all samples in imaging_request_1 in request_name: zimmerman_01 then do:
# 
# clearmap_pipeline.sh zimmerman_01 imaging_request_1
# 
# If you want to run for just sample zimmerman_01-002 then do:
# 
# clearmap_pipeline.sh zimmerman_01 imaging_request_1 zimmerman_01-002
###
###
### NOTE: This script will only work on one imaging request at a time. 
###
### NOTE2: Set the username, output_rootpath, clearmap_params_file and atlas below
### You can make a new clearmap params pickle file of custom parameters using the script:
### spock-clearmap/clearmap_parameters.py
username=$1
output_rootpath="/jukebox/wang/sanjeev/cm_output"
clearmap_params_file='/jukebox/wang/sanjeev/cm_output/cell_detection_parameter.p'
atlas='Allen' # 'Princeton' or 'Allen'
# Makes sure that you supplied at least two command line arguments 
if [ "$#" -lt 2 ]
then
	echo "Error: Incorrect number of command line arguments"
	echo ""
	echo "Usage: clearmap_pipeline.sh request_name imaging_request"
	echo "e.g.: clearmap_pipeline.sh zimmerman_01 imaging_request_1"
	exit 1
fi

request_name=$2
imaging_request=$3 
request_dir="/jukebox/LightSheetData/lightserv/${username}/${request_name}"

# Check to see if sample_name argument was provided 
if [ "$#" -eq 4 ]
then
	sample_name=$4
	echo "Running script for Request name: $2, Imaging request: $3, Sample name: $4"
	declare -a arr=("/jukebox/LightSheetData/lightserv/${username}/${request_name}/${sample_name}")
else
	echo "Running script for Request name: $2, Imaging request: $3, all samples"
	declare -a arr=(${request_dir}/*)
fi
echo ""
# echo $arr
echo "Sample directories that will be run are:"
for d in "${arr[@]}"
do
	echo $d
done
echo ""
# Loop through all sample dirs in the request folder and run the full pipeline for each 

blocks_per_job=1 # number of cell blocks to process in a single array job. Currently only works with 1 right now
# due to slurm error. TODO - allow parallel processing in each array job 

# for sample_dir in "${request_dir}"/*
for sample_dir in "${arr[@]}"
do
	echo "Working on sample: $sample_dir"
	# Step 0: Link over raw 488 and 642 files to destination directory.
	# Also creates a JSON file (the blockfile) which just contains the number of blocks to do cell detection
	echo "Step 0: Synchronous preprocessing"
	module load anacondapy/2020.11
	conda activate ClearMap
	python spock-clearmap/clearmap_preprocessing.py ${sample_dir} ${imaging_request} ${output_rootpath} 2>&1 | tee
	conda deactivate 
	module unload anacondapy/2020.11
	sample_name=$(basename ${sample_dir})
	blockfile="${output_rootpath}/${request_name}/${sample_name}/${imaging_request}/rawdata/resolution_3.6x/block_processing_info.json"		
	
	# Read JSON file to figure out how many blocks there are to process
	if [ -a $blockfile ]
	then
		# echo $blockfile
		n_blocks=`cat ${blockfile} | cut -d " " -f 2 | cut -d "}" -f 1`
		echo "Have ${n_blocks} blocks for clearmap to process"
	else
		echo "Block file not found. Skipping this sample"
		echo ""
		continue
	fi
	max_array_job=`echo "(${n_blocks}+${blocks_per_job}-1)/${blocks_per_job}-1" | bc` # bash equivalent of ceil() 
	# echo "Max array job for cell detection: ${max_array_job}"
	echo "Submitting batch jobs:" 

	## Create stitched.npy memmap volume file 
	OUT1=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/create_memmap_vol.sh)
	echo "Step 1: Memmmap volume step"
	echo $OUT1

	# # ## Run the individual blocks in array jobs
	OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } \
	--exclude=./bad_nodenames.txt \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},blocks_per_job=${blocks_per_job},\
output_rootpath=${output_rootpath},clearmap_params_file=${clearmap_params_file} \
	--array=0-${max_array_job} spock-clearmap/slurm_scripts/cell_detect.sh)
	echo "Step 2: Cell detection on blocks:"
	echo $OUT2
	
	# ## Merge the blocks
	OUT3=$(sbatch --parsable --dependency=afterok:${OUT2##* } \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},\
clearmap_params_file=${clearmap_params_file} \
	spock-clearmap/slurm_scripts/merge_blocks.sh)
	echo "Step 3: Merge cell detection blocks:"
	echo $OUT3
	
done