#!/bin/env bash

### How to run:
# clearmap_pipeline.sh ${request_name}  

# echo "Request name: $1"

# Loop through all sample dirs in the request folder and run the full pipeline for each 
request_name=$1
request_dir="/jukebox/LightSheetData/lightserv/cz15/$1"
output_rootpath="/jukebox/wang/ahoag/for_cz/clearmap2_test_output"
declare -a arr=("/jukebox/wang/ahoag/for_cz/clearmap2_test_output/zimmerman_01/zimmerman_01-002")
# for sample_dir in "${request_dir}"/*
for sample_dir in "${arr[@]}"
do
	echo $sample_dir
	# # First link over raw 488 and 642 files to destination directory
	# echo "Running preprocessing step synchronously"
	# python spock-clearmap/cz_clearmap_preprocessing.py ${sample_dir} 2>&1 | tee

	# sample_name=$(basename ${sample_dir})
	# blockfile="${output_rootpath}/$1/${sample_name}/imaging_request_1/rawdata/resolution_3.6x/block_processing_info.json"		
	
	# # Read JSON file to figure out how many blocks there are to process
	# if [ -a $blockfile ]
	# then
	# 	echo $blockfile
	# 	n_blocks=`cat ${blockfile} | cut -d " " -f 2 | cut -d "}" -f 1`
	# 	echo ${n_blocks}
	# else
	# 	echo "Block file not found"
	# 	continue
	# fi
	# max_array_job=`echo "$n_blocks-1" | bc`
	# echo "Max array job: $max_array_job"
	### Create stitched.npy memmap volume file 
	# OUT1=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir} \
	# spock-clearmap/slurm_scripts/create_memmap_vol_test.sh)
	# echo $OUT1

	### Run the individual blocks in array jobs
	# OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --exclude=./bad_nodenames.txt \
	# --export=ALL,sample_dir=${sample_dir} --array=40-41 spock-clearmap/slurm_scripts/cell_detect_test.sh)
	# OUT2=$(sbatch --parsable --exclude=./bad_nodenames.txt \
	# --export=ALL,sample_dir=${sample_dir} --array=40-41 spock-clearmap/slurm_scripts/cell_detect_test.sh)
	# echo $OUT2
	
	### Merge the blocks
	# OUT3=$(sbatch --parsable --dependency=afterok:${OUT2##* } \
	# --export=ALL,sample_dir=${sample_dir} spock-clearmap/slurm_scripts/merge_blocks_test.sh)
	OUT3=$(sbatch --parsable \
	--export=ALL,sample_dir=${sample_dir} spock-clearmap/slurm_scripts/merge_blocks_test.sh)
done



# #Cell Detection
# OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --export=ALL,request_name=$2 spock-clearmap/slurm_scripts/main.sh)
# echo $OUT2

# #Downsizing
# OUT3=$(sbatch --parsable --dependency --export=ALL,corrected_dir=$1,downsized_dir=$3 downsizing/spim_downsize.sh)
# echo $OUT3

# #Inverse registration
# OUT4=$(sbatch --parsable --dependency=afterok:${OUT3##* } --export=ALL,src=$3,reg=$4,cell=$5 registration/slurm_scripts/spim_inverse_register.sh)
# echo $OUT4

# #Register cells and make CSV data frame of counts in each region
# OUT5=$(sbatch --parsable --dependency=afterok:${OUT2##* }:${OUT4##* } --export=ALL,request_name=$2 spock-clearmap/slurm_scripts/postprocessing.sh)
# echo $OUT5

