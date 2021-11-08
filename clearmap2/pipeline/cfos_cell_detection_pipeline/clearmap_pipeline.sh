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
username=cz15
output_rootpath="/jukebox/wang/ahoag/for_cz/clearmap2_test_output"
clearmap_params_file='/jukebox/witten/Chris/data/clearmap2/utilities/cell_detection_parameter.p'
atlas='Princeton' # 'Princeton' or 'Allen'
# Makes sure that you supplied at least two command line arguments 
if [ "$#" -lt 2 ]
then
	echo "Error: Incorrect number of command line arguments"
	echo ""
	echo "Usage: clearmap_pipeline.sh request_name imaging_request"
	echo "e.g.: clearmap_pipeline.sh zimmerman_01 imaging_request_1"
	exit 1
fi

request_name=$1
imaging_request=$2 
request_dir="/jukebox/LightSheetData/lightserv/${username}/$1"

# Check to see if sample_name argument was provided 
if [ "$#" -eq 3 ]
then
	sample_name=$3
	echo "Running script for Request name: $1, Imaging request: $2, Sample name: $3"
	declare -a arr=("/jukebox/LightSheetData/lightserv/${username}/${request_name}/${sample_name}")
else
	echo "Running script for Request name: $1, Imaging request: $2, all samples"
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
	blockfile="${output_rootpath}/$1/${sample_name}/${imaging_request}/rawdata/resolution_3.6x/block_processing_info.json"		
	
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
	
	## Diagnostic plot for corrected planes
	OUT0_dg=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/diagnostic_corrected_planes.sh)
	echo "Step 0-diag: Diagnostic plots of corrected planes:"
	echo $OUT0_dg

	## Create stitched.npy memmap volume file 
	OUT1=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/create_memmap_vol.sh)
	echo "Step 1: Memmmap volume step"
	echo $OUT1

	## Diagnostic plot for memmap vol
	OUT1_dg=$(sbatch --parsable --dependency=afterok:${OUT1##* } \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/diagnostic_memmap_vol.sh)
	echo "Step 1-diag: Diagnostic plots of memmap volume:"
	echo $OUT1_dg

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

	# Diagnostic plot for merged cells
	OUT3_dg=$(sbatch --parsable --dependency=afterok:${OUT3##* } \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/diagnostic_merged_cells.sh)
	echo "Step 3-diag: Diagnostic plots of merged cells:"
	echo $OUT3_dg

	## Downsizing, both channels one per array job, can start without dependency
	OUT4=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},\
atlas=${atlas} --array=0-1 downsizing/spim_downsize.sh)
	echo "Step 4: Downsizing:"
	echo $OUT4

	## Diagnostic plot for downsized sagittal planes
	OUT4_dg=$(sbatch --parsable --dependency=afterok:${OUT4##* } \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath} \
	downsizing/diagnostic_downsize.sh)
	echo "Step 4-diag: Diagnostic plots of downsized planes:"
	echo $OUT4_dg

	## Inverse registration, both transformations one per array job
	OUT5=$(sbatch --parsable --dependency=afterok:${OUT4##* } --array=0-1 \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},\
atlas=${atlas} registration/slurm_scripts/spim_inverse_register.sh)
	echo "Step 5: Inverse registration:"
	echo $OUT5

	## Register cells to atlas space and make CSV data frame of counts in each region
	# Dependent on merge block step and inverse registration step
	OUT6=$(sbatch --parsable --dependency=afterok:${OUT3##* }:${OUT5##* } \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath},atlas=${atlas} \
	spock-clearmap/slurm_scripts/postprocessing.sh)
	echo "Step 6: Register cells to atlas and create brain region count dataframe"
	echo $OUT6

	## Diagnostic plot for final data products (registered cells and brain region count dataframe)
	OUT6_dg=$(sbatch --parsable --dependency=afterok:${OUT6##* } \
	--export=ALL,sample_dir=${sample_dir},\
imaging_request=${imaging_request},output_rootpath=${output_rootpath} \
	spock-clearmap/diagnostic_postprocessing.sh)
	echo "Step 6-diag: Diagnostic plots of registered cells and brain region count dataframe:"
	echo $OUT6_dg
	echo ""
done



