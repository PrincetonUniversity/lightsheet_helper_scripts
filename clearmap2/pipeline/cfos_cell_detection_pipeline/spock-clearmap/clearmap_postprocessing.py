# General imports
import os,sys,pickle
import numpy as np
import pandas as pd
import json
import tifffile
from concurrent.futures import ProcessPoolExecutor
from brain_atlas_toolkit import graph_tools

# ClearMap2 imports
sys.path.append('/jukebox/braininit/lightsheet/ClearMap2')

# os.chmod("/jukebox/witten/Chris/python/ClearMap2-master/ClearMap/External/elastix/build/bin/transformix",0o777)
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ImageProcessing.Experts.Cells as cells
import ClearMap.Settings as settings
import ClearMap.Alignment.Resampling as res
import ClearMap.Alignment.Elastix as elx

from functools import partial




def get_count_and_volume(region_idx,segment_name_dict,eroded_atlas_vol,segment_list,):
	""" 
	---PURPOSE---
	Given a list of segments, get the number of counts in each brain region.
	---INPUT---
	region_idx        - 1D array containing the brain region ID in which each cell was detected 
	segment_name_dict - dictionary mapping ids:segment names
	eroded_atlas_vol  - the eroded atlas volume
	segment_list      - list of integer ids
	---OUTPUT---
	count_dict   - dictionary where keys are segment ids and values are counts
	volume_dict  - dictionary where keys are segment ids and values are volume of that region in the eroded atlas
	"""
	count_dict = {}
	volume_dict = {}
	print(f"parallel processing {segment_list}")
	for atlas_segment_id in segment_list:
		segment_name = segment_name_dict[int(atlas_segment_id)]
		if int(atlas_segment_id) not in atlas_segments:
			count_dict[segment_name] = 0
			volume_dict[segment_name] = 0
		else:
			count_dict[segment_name] = np.count_nonzero(region_idx==int(atlas_segment_id))
			volume_dict[segment_name] = np.sum(eroded_atlas_vol == int(atlas_segment_id)) * np.power(0.02,3)
	return count_dict,volume_dict

if __name__ == "__main__":
	n_cores = os.cpu_count()
	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")
	atlas = sys.argv[4].strip().rstrip("/")

	if atlas == 'Princeton':
		# Princeton Mouse Atlas
		eroded_atlas_file = '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_16bit_hierarch_labels_60um_edge_80um_vent_erosion.tif'
		segment_props_file = '/jukebox/LightSheetTransfer/atlas/PMA_16bit_hierarch_labels_segment_properties_info'
		ontology_json_file = '/jukebox/LightSheetTransfer/atlas/PMA_ontology.json'
	elif atlas == 'Allen':
		# Allen Mouse Brain Atlas
		eroded_atlas_file = '/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_60um_edge_80um_vent_erosion.tif'
		segment_props_file = '/jukebox/LightSheetTransfer/atlas/allen_atlas/allenatlas_2017_16bit_hierarch_labels_segment_properties_info'
		ontology_json_file = '/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json'
	else:
		sys.exit(f"Atlas provided: {atlas} is not accepted. Must be one of ['Princeton','Allen']")


	request_name,sample_name = sample_dir.split('/')[-2:]
	workspace_dir =  os.path.join(output_rootpath,
		request_name,sample_name,imaging_request,
		"rawdata/resolution_3.6x")
	
	#  Set paths to Elastix transformation files and downsized files
	elastix_inverse_dir = os.path.join(workspace_dir,'elastix_inverse_transform')
	
	ch488_downsized_dir = os.path.join(workspace_dir,'Ex_488_Em_0_downsized')
	ch642_downsized_dir = os.path.join(workspace_dir,'Ex_642_Em_2_downsized')
	ch488_downsized_file = os.path.join(ch488_downsized_dir,'downsized_for_atlas_ch488.tif')
	ch642_downsized_file = os.path.join(ch642_downsized_dir,'downsized_for_atlas_ch642.tif')

	# Initialize ClearMap2 workspace object
	ws = wsp.Workspace('CellMap',directory=workspace_dir)
	ws.debug = False
	print()
	ws.info()

	# Load detected cells from cells_raw.npy
	raw_source = ws.source('cells', postfix='raw')
	size_intensity = np.hstack([raw_source[c][:,None] for c in ['size','background']])
	coordinates_raw = np.hstack([raw_source[c][:,None] for c in 'xyz'])
	
	# Swap axes to go from horizontal to sagittal orientation
	coordinates_raw_swapped_axes = np.zeros_like(coordinates_raw)
	coordinates_raw_swapped_axes[:,0] = coordinates_raw[:,2]
	coordinates_raw_swapped_axes[:,1] = coordinates_raw[:,1]
	coordinates_raw_swapped_axes[:,2] = coordinates_raw[:,0]

	# Transform cell coordinates from raw 642-space to downsampled 642-space
	coordinates_resampled = res.resample_points(
					coordinates_raw_swapped_axes, sink=None, orientation=None,
					source_shape=io.shape(ws.filename('stitched'))[::-1],
					sink_shape=io.shape(ch642_downsized_file))

	# Transform cell coordinates from 642-space to 488-space
	coordinates_aligned_to_488 = elx.transform_points(
					coordinates_resampled, sink=None,
					transform_directory=os.path.join(elastix_inverse_dir,
						'488_to_642'),
					temp_file='/tmp/elastix_input_pipeline.bin',
					result_directory='/tmp/elastix_output_pipeline')

	# Change permissions back since transformix alters permissions when it is run
	os.chmod(os.path.join(elastix_inverse_dir,
		'488_to_642','TransformParameters.0.txt'),0o777)
	os.chmod(os.path.join(elastix_inverse_dir,
		'488_to_642','TransformParameters.1.txt'),0o777)

	# Transform cell coordinates from 488-space to atlas-space
	coordinates_aligned_to_atlas = elx.transform_points(
					coordinates_aligned_to_488, sink=None,
					transform_directory=elastix_inverse_dir,
					binary=True, indices=False,
					temp_file='/tmp/elastix_input_testcz15.bin',
					result_directory='/tmp/elastix_output_cztest')

	# Change permissions back since transformix alters permissions when it is run
	os.chmod(os.path.join(elastix_inverse_dir,'TransformParameters.0.txt'),0o777)
	os.chmod(os.path.join(elastix_inverse_dir,'TransformParameters.1.txt'),0o777)

	# Load atlas files
	eroded_atlas_vol = np.array(tifffile.imread(eroded_atlas_file)).astype('uint16')
	atlas_segments = np.unique(eroded_atlas_vol)
	atlas_segments = np.array([x for x in atlas_segments if x!=0])
	
	with open(segment_props_file,'r') as infile:
		segment_props_dict = json.load(infile)
	
	with open(ontology_json_file,'r') as infile:
		ontology_dict = json.load(infile)

	# Record the brain region ID where a cell is detected, 0 if not in a region
	cell_regions = np.empty([len(coordinates_aligned_to_atlas), 1], dtype=int)
	xyz = np.asarray([(int(X[0]), int(X[1]), int(X[2])) for X in coordinates_aligned_to_atlas])
	for idx, val in enumerate(xyz):
		try:
			ID = eroded_atlas_vol[val[2],val[1],val[0]]
			cell_regions[idx] = ID
		except Exception as e:
			cell_regions[idx] = 0
			pass
	
	# Add brain region ID to transformed cell array 
	cells_to_save = np.hstack((coordinates_aligned_to_atlas,size_intensity,cell_regions))
	header = ['x','y','z','size','intensity','region']
	dtypes = [int, int, int, int, float, int]
	dt = {'names' : header, 'formats' : dtypes}
	output_array = np.zeros(len(cells_to_save), dtype=dt)
	for i,h in enumerate(header):
		output_array[h] = cells_to_save[:,i]
	# Remove cells that are outside the atlas
	output_array = np.delete(output_array,np.argwhere(cell_regions==0))

	# Save registered cells to cells_transformed_to_atlas.npy
	savename = ws.filename('cells',postfix='transformed_to_atlas')
	io.write(savename,output_array)
	print(f'Saving registered cell detection results to: {savename}')
	print()

	# Filter cells and save to cells_transformed_to_atlas_filtered.npy
	thresholds_file = '/jukebox/witten/Chris/data/clearmap2/utilities/cell_detection_filter.p'
	with open(thresholds_file,'rb') as f:
		thresholds_dict = pickle.load(f)
	minsize_thresh = str(thresholds_dict['size'][0])
	postfix_filtered = f"transformed_to_atlas_filtered_{minsize_thresh}px"
	cells_filtered = cells.filter_cells(
					 source = ws.filename('cells', postfix='transformed_to_atlas'),
					 sink = ws.filename('cells', postfix=postfix_filtered),
					 thresholds=thresholds_dict)
	print(f"Saving filtered cell detection results to: {ws.filename('cells',postfix=postfix_filtered)}")
	print()

	# Get cell counts for each brain region
	region_idx = ws.source('cells', postfix=postfix_filtered)['region']
	ids = segment_props_dict['inline']['ids']
	segment_names = segment_props_dict['inline']['properties'][0]['values']
	segment_name_dict = {int(ids[ii]):segment_names[ii].split(':')[1].strip() for ii in range(len(ids))}
	count_dict_names = {}
	volume_dict_names = {}
	# Loop over all ids in the atlas JSON, not just the ones in the volume
	print("Determining counts and volumes for each region")
	count_dict_names = {}
	volume_dict_names = {}
	chunk_size = 25 # Each core runs get_count_and_volume() on this many different regions
	chunked_segment_lists = [ids[i:i+chunk_size] for i in range(0,len(ids),chunk_size)]
	get_count_and_volume_par = partial(get_count_and_volume,
		region_idx,segment_name_dict,eroded_atlas_vol)
	with ProcessPoolExecutor(max_workers=n_cores) as executor:
		for count_dict_i,volume_dict_i in executor.map(get_count_and_volume_par,chunked_segment_lists):
			try:
				for key in count_dict_i:
					count_dict_names[key] = count_dict_i[key]
				for key in volume_dict_i:
					volume_dict_names[key] = volume_dict_i[key]
			except Exception as exc:
				print(f'generated an exception: {exc}')
	sys.stdout.flush()
	print()
	print("Correcting counts")
	# Correct counts by adding sum of counts in progeny regions
	ontology_graph = graph_tools.Graph(ontology_dict)
	corrected_count_dict = {}
	corrected_density_dict = {}
	for region in count_dict_names.keys():
		counts_region = count_dict_names[region]
		volume_region = volume_dict_names[region]
		progeny = ontology_graph.get_progeny(region)
		if progeny != []:
			for prog in progeny:
				try:
					counts_region += count_dict_names[prog]
					volume_region += volume_dict_names[prog]
				except KeyError:
					continue
		corrected_count_dict[region] = counts_region
		try:
			corrected_density_dict[region] = counts_region/volume_region
		except ZeroDivisionError:
			corrected_density_dict[region] = float('NaN')
	print("Corrected counts. Now saving out final CSV file")
	sys.stdout.flush()
	# Save region cell counts to region_cell_counts.csv
	df = pd.DataFrame([corrected_count_dict,corrected_density_dict])
	basename_csv = f"region_cell_counts_filtered_{minsize_thresh}px.csv"
	savename_csv = os.path.join(workspace_dir,basename_csv)
	df.to_csv(savename_csv,index=False)
	print(f'Saving cell detection results broken down by brain region to: {savename_csv}')
	print()
