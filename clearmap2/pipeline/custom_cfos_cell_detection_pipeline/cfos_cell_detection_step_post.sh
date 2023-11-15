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
output_rootpath="/jukebox/wang/sanjeev/cm_output/"
clearmap_params_file='/jukebox/wang/sanjeev/cm_output/cell_detection_parameter.p'
atlas='Allen' # 'Princeton' or 'Allen'

request_dir=$1
request_name=$2
microscope=$3

# Check to see if sample_name argument was provided 
echo "Running script for Request name: $1, all samples"
declare -a arr=(${request_dir}/*)

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
	python spock-clearmap/clearmap_preprocessing.py ${sample_dir} ${request_name} ${output_rootpath} ${microscope} 2>&1 | tee
	conda deactivate 
	module unload anacondapy/2020.11
	sample_name=$(basename ${sample_dir})
	blockfile="${output_rootpath}/$2/${sample_name}/block_processing_info.json"	

	echo $blockfile	
	
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

	OUT1=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},\
request_name=${request_name},output_rootpath=${output_rootpath},atlas=${atlas} \
	spock-clearmap/slurm_scripts/postprocessing.sh)
	echo "Step 6: Register cells to atlas and create brain region count dataframe"
	echo $OUT1

done