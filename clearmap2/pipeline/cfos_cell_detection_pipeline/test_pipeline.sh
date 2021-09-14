#!/bin/env bash

### How to run:
# clearmap_pipeline.sh ${request_name}  
target="/jukebox/LightSheetData/lightserv/cz15/$1"
output_rootpath="/jukebox/wang/ahoag/for_cz/clearmap2_test_output"
for sample_dir in "${target}"/*
do
	sample_name=$(basename ${sample_dir})
	blockfile="${output_rootpath}/$1/${sample_name}/imaging_request_1/rawdata/resolution_3.6x/block_processing_info.json"
	if [ -a $blockfile ]
	then
		echo $blockfile
		n_blocks=`cat ${blockfile} | cut -d " " -f 2 | cut -d "}" -f 1`
		echo ${n_blocks}
	else
		echo "No block file found, skipping rest of pipeline for this sample "
		continue
	fi
done
