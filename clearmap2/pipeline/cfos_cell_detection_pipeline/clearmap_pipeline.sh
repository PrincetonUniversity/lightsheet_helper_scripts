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
	echo "Step 0: Synchronous preprocessing"
	module load anacondapy/2020.11
	conda activate ClearMap
	python spock-clearmap/cz_clearmap_preprocessing.py ${sample_dir} ${output_rootpath} 2>&1 | tee
	conda deactivate 
	module unload anacondapy/2020.11
	sample_name=$(basename ${sample_dir})
	blockfile="${output_rootpath}/$1/${sample_name}/imaging_request_1/rawdata/resolution_3.6x/block_processing_info.json"		
	
	# Read JSON file to figure out how many blocks there are to process
	if [ -a $blockfile ]
	then
		# echo $blockfile
		n_blocks=`cat ${blockfile} | cut -d " " -f 2 | cut -d "}" -f 1`
		# echo ${n_blocks}
	else
		echo "Block file not found"
		exit
	fi
	max_array_job=`echo "(${n_blocks}+${blocks_per_job}-1)/${blocks_per_job}-1" | bc` # bash equivalent of ceil() 
	echo "Max array job for cell detection: ${max_array_job}"
	echo "Submitting batch jobs:" 

	## Diagnostic plot for corrected planes
	OUT0_dg=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/diagnostic_corrected_planes.sh)
	echo "Step 0-diag: Diagnostic plots of corrected planes:"
	echo $OUT0_dg

	## Create stitched.npy memmap volume file 
	OUT1=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/create_memmap_vol.sh)
	echo "Step 1: Memmmap volume step"
	echo $OUT1

	## Diagnostic plot for corrected planes
	OUT1_dg=$(sbatch --parsable --dependency=afterok:${OUT1##* } \
	--export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/diagnostic_memmap_vol.sh)
	echo "Step 1-diag: Diagnostic plots of memmap volume:"
	echo $OUT1_dg

	# # ## Run the individual blocks in array jobs
	OUT2=$(sbatch --parsable --dependency=afterok:${OUT1##* } --exclude=./bad_nodenames.txt \
	--export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath},blocks_per_job=${blocks_per_job} \
	--array=0-${max_array_job} spock-clearmap/slurm_scripts/cell_detect.sh)
	echo "Step 2: Cell detection on blocks:"
	echo $OUT2
	
	# ## Merge the blocks
	OUT3=$(sbatch --parsable --dependency=afterok:${OUT2##* } \
	--export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/merge_blocks.sh)
	echo "Step 3: Merge cell detection blocks:"
	echo $OUT3

	## Diagnostic plot for merged cells
	OUT3_dg=$(sbatch --parsable --dependency=afterok:${OUT3##* } \
	--export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/diagnostic_merged_cells.sh)
	echo "Step 3-diag: Diagnostic plots of merged cells:"
	echo $OUT3_dg

	## Downsizing, both channels one per array job, can start without dependency
	OUT4=$(sbatch --parsable --export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	--array=0-1 downsizing/spim_downsize.sh)
	echo "Step 4: Downsizing:"
	echo $OUT4

	## Diagnostic plot for downsized sagittal planes
	OUT4_dg=$(sbatch --parsable --dependency=afterok:${OUT4##* } \
	--export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	downsizing/diagnostic_downsize.sh)
	echo "Step 4-diag: Diagnostic plots of downsized planes:"
	echo $OUT4_dg

	## Inverse registration, both transformations one per array job
	OUT5=$(sbatch --parsable --dependency=afterok:${OUT4##* } --array=0-1 \
	--export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	registration/slurm_scripts/spim_inverse_register.sh)
	echo "Step 5: Inverse registration:"
	echo $OUT5

	## Register cells to atlas space and make CSV data frame of counts in each region
	# Dependent on merge block step and inverse registration step
	OUT6=$(sbatch --parsable --dependency=afterok:${OUT3##* }:${OUT5##* } \
	--export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	spock-clearmap/slurm_scripts/postprocessing.sh)
	echo "Step 6: Register cells to atlas and create brain region count dataframe"
	echo $OUT6

	## Diagnostic plot for final data products (registered cells and brain region count dataframe)
	OUT6_dg=$(sbatch --parsable --dependency=afterok:${OUT6##* } \
	--export=ALL,sample_dir=${sample_dir},output_rootpath=${output_rootpath} \
	spock-clearmap/diagnostic_postprocessing.sh)
	echo "Step 6-diag: Diagnostic plots of registered cells and brain region count dataframe:"
	echo $OUT6_dg

done



