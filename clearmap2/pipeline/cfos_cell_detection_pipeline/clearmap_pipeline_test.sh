#!/bin/env bash

### How to run:
# clearmap_pipeline.sh ${request_name}  

# echo "Request name: $1"

# Loop through all sample dirs in the request folder and run the full pipeline for each 
request_name=$1
request_dir="/jukebox/LightSheetData/lightserv/cz15/$1"
output_rootpath="/jukebox/wang/ahoag/for_cz/clearmap2_test_output"
blocks_per_job=2 # number of cell blocks to process in a single array job
declare -a arr=("/jukebox/LightSheetData/lightserv/cz15/zimmerman_01/zimmerman_01-002")
# for sample_dir in "${request_dir}"/*
for sample_dir in "${arr[@]}"
do
	echo $sample_dir
	# First link over raw 488 and 642 files to destination directory
	echo "Running preprocessing step synchronously"
	module load anacondapy/2020.11
	conda activate ClearMap
	python spock-clearmap/cz_clearmap_preprocessing.py ${sample_dir} 2>&1 | tee
	conda deactivate 
	module unload anacondapy/2020.11
	sample_name=$(basename ${sample_dir})
	blockfile="${output_rootpath}/$1/${sample_name}/imaging_request_1/rawdata/resolution_3.6x/block_processing_info.json"		
	
	# Read JSON file to figure out how many blocks there are to process
	if [ -a $blockfile ]
	then
		echo $blockfile
		n_blocks=`cat ${blockfile} | cut -d " " -f 2 | cut -d "}" -f 1`
		echo ${n_blocks}
	else
		echo "Block file not found"
		exit
	fi
	max_array_job=`echo "(${n_blocks}+${blocks_per_job}-1)/${blocks_per_job}-1" | bc` # bash equivalent of ceil() 
	echo "Max array job for cell detection: ${max_array_job}"
	echo "Submitting batch jobs:" 
	## Create stitched.npy memmap volume file 
	OUT1=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir} \
	spock-clearmap/slurm_scripts/create_memmap_vol_test.sh)
	echo $OUT1

	# # ## Run the individual blocks in array jobs
	OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --exclude=./bad_nodenames.txt \
	--export=ALL,sample_dir=${sample_dir},blocks_per_job=${blocks_per_job} --array=0-${max_array_job} \
	spock-clearmap/slurm_scripts/cell_detect_test.sh)
	echo $OUT2
	
	# ## Merge the blocks
	OUT3=$(sbatch --parsable --dependency=afterok:${OUT2##* } \
	--export=ALL,sample_dir=${sample_dir} spock-clearmap/slurm_scripts/merge_blocks_test.sh)
	echo $OUT3

	## Downsizing, both channels one per array job, can start without dependency
	OUT4=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir} --array=0-1 downsizing/spim_downsize_test.sh)
	echo $OUT4

	## Inverse registration, both transformations one per array job
	OUT5=$(sbatch --parsable --dependency=afterok:${OUT4##* } --array=0-1 \
	--export=ALL,sample_dir=${sample_dir} registration/slurm_scripts/spim_inverse_register_test.sh)
	echo $OUT5

	## Register cells to atlas space and make CSV data frame of counts in each region
	# Dependent on merge block step and inverse registration step
	OUT6=$(sbatch --parsable --dependency=afterok:${OUT3##* }:${OUT5##* } \
	--export=ALL,sample_dir=${sample_dir} spock-clearmap/slurm_scripts/postprocessing_test.sh)
	echo $OUT6

done



