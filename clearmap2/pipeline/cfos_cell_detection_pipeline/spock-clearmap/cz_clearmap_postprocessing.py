# General imports
import os,sys,pickle
import numpy as np
import pandas as pd
import json
import tifffile
from concurrent.futures import ProcessPoolExecutor
from brain_atlas_toolkit import graph_tools

# ClearMap2 imports
sys.path.append('/jukebox/witten/Chris/python/ClearMap2-master')
os.chmod("/jukebox/witten/Chris/python/ClearMap2-master/ClearMap/External/elastix/build/bin/transformix",0o777)
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ImageProcessing.Experts.Cells as cells
import ClearMap.Settings as settings
import ClearMap.Alignment.Resampling as res
import ClearMap.Alignment.Elastix as elx
import ClearMap.Utils.HierarchicalDict as hdict

# Select dataset to analyze
fpath = os.path.join(sys.argv[1],'rawdata/resolution_3.6x')

#  Set paths to Elastix transformation files
data_dir = os.path.join('/jukebox/LightSheetData/lightserv/cz15',fpath)
elastix_inverse_dir = os.path.join(data_dir,'elastix_inverse_transform')
ch488_dir = os.path.join(data_dir,'Ex_488_Em_0_downsized')
ch642_dir = os.path.join(data_dir,'Ex_642_Em_2_downsized')
ch488_downsized_file = os.path.join(ch488_dir,'downsized_for_atlas.tif')
ch642_downsized_file = os.path.join(ch642_dir,'downsized_for_atlas.tif')

# Initialize ClearMap2 workspace object
directory = os.path.join('/jukebox/witten/Chris/data/clearmap2',fpath)
expression_raw = '/Ex_642_Em_2_corrected/Z<Z,4>.tif'
expression_auto = '/Ex_488_Em_0_corrected/Z<Z,4>.tif'
ws = wsp.Workspace('CellMap',directory=directory)
ws.update(raw=expression_raw)
ws.update(autofluorescence=expression_auto)
resources_directory=settings.resources_path
ws.debug = False
print()
ws.info()

# Load detected cells from cells_raw.npy
raw_source = ws.source('cells', postfix='raw')
size_intensity = np.hstack([raw_source[c][:,None] for c in ['size','background']])
coordinates_raw = np.hstack([raw_source[c][:,None] for c in 'xyz'])
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
                  transform_directory=os.path.join(elastix_inverse_dir,'Ex_488_Em_0_downsized_to_Ex_642_Em_2_downsized'),
                  temp_file='/tmp/elastix_input_testcz15.bin',
                  result_directory='/tmp/elastix_output_cztest')
os.chmod(os.path.join(elastix_inverse_dir,'Ex_488_Em_0_downsized_to_Ex_642_Em_2_downsized','TransformParameters.0.txt'),0o777)
os.chmod(os.path.join(elastix_inverse_dir,'Ex_488_Em_0_downsized_to_Ex_642_Em_2_downsized','TransformParameters.1.txt'),0o777)

# Transform cell coordinates from 488-space to atlas-space
coordinates_aligned_to_atlas = elx.transform_points(
                  coordinates_aligned_to_488, sink=None,
                  transform_directory=elastix_inverse_dir,
                  binary=True, indices=False,
                  temp_file='/tmp/elastix_input_testcz15.bin',
                  result_directory='/tmp/elastix_output_cztest')
os.chmod(os.path.join(elastix_inverse_dir,'TransformParameters.0.txt'),0o777)
os.chmod(os.path.join(elastix_inverse_dir,'TransformParameters.1.txt'),0o777)

# Load atlas files
eroded_atlas_file = '/jukebox/witten/Chris/data/clearmap2/utilities/princeton-atlas/annotation_sagittal_atlas_20um_16bit_hierarch_labels_60um_edge_80um_vent_erosion.tif'
eroded_atlas_vol = np.array(tifffile.imread(eroded_atlas_file)).astype('uint16')
atlas_segments = np.unique(eroded_atlas_vol)
atlas_segments = np.array([x for x in atlas_segments if x!=0])
segment_props_file = '/jukebox/witten/Chris/data/clearmap2/utilities/princeton-atlas/pma_segment_properties_info'
with open(segment_props_file,'r') as infile:
    segment_props_dict = json.load(infile)
ontology_json_file = '/jukebox/witten/Chris/data/clearmap2/utilities/princeton-atlas/PMA_ontology.json'
with open(ontology_json_file,'r') as infile:
    ontology_dict = json.load(infile)

# Remove cells that are outside the atlas
cell_regions = np.empty([len(coordinates_aligned_to_atlas), 1], dtype=int)
xyz = np.asarray([(int(X[0]), int(X[1]), int(X[2])) for X in coordinates_aligned_to_atlas])
for idx, val in enumerate(xyz):
    try:
        ID = eroded_atlas_vol[val[2],val[1],val[0]]
        cell_regions[idx] = ID
    except Exception as e:
        cell_regions[idx] = 0
        pass
cells_to_save = np.hstack((coordinates_aligned_to_atlas,size_intensity,cell_regions))
header = ['x','y','z','size','intensity','region']
dtypes = [int, int, int, int, float, int]
dt = {'names' : header, 'formats' : dtypes}
output_array = np.zeros(len(cells_to_save), dtype=dt)
for i,h in enumerate(header):
    output_array[h] = cells_to_save[:,i]
output_array = np.delete(output_array,np.argwhere(cell_regions==0))

# Save registered cells to cells_transformed_to_atlas.npy
savename = ws.filename('cells',postfix='transformed_to_atlas')
io.write(savename,output_array)
print('Saving registered cell detection results to:\n/jukebox/witten/Chris/data/clearmap2/' + fpath + '/cells_transformed_to_atlas.npy')
print()

# Filter cells and save to cells_transformed_to_atlas_filtered.npy
fname = '/jukebox/witten/Chris/data/clearmap2/utilities/cell_detection_filter.p'
with open(fname,'rb') as f:
    cell_detection_filter = pickle.load(f)
print('path     : ' + fname)
hdict.pprint(cell_detection_filter)
print()
thresholds = {'x':cell_detection_filter['x'],'y':cell_detection_filter['y'],'z':cell_detection_filter['z'],'size':cell_detection_filter['size'],'intensity':cell_detection_filter['intensity'],'region':cell_detection_filter['region']}
fname = 'transformed_to_atlas_filtered_' + str(cell_detection_filter['size'][0]) + 'px'
cells_filtered = cells.filter_cells(
                   source = ws.filename('cells', postfix='transformed_to_atlas'),
                   sink = ws.filename('cells', postfix=fname),
                   thresholds=thresholds)
print('Saving filtered cell detection results to:\n/jukebox/witten/Chris/data/clearmap2/' + fpath + '/cells_' + fname + '.npy')
print()

# Get cell counts for each brain region
region_idx = ws.source('cells', postfix=fname)['region']
ids = segment_props_dict['inline']['ids']
segment_names = segment_props_dict['inline']['properties'][0]['values']
segment_name_dict = {int(ids[ii]):segment_names[ii].split(':')[1].strip() for ii in range(len(ids))}
count_dict_names = {}
volume_dict_names = {}
for atlas_segment_id in ids:
    segment_name = segment_name_dict[int(atlas_segment_id)]
    if int(atlas_segment_id) not in atlas_segments:
        count_dict_names[segment_name] = 0
        volume_dict_names[segment_name] = 0
    else:
        count_dict_names[segment_name] = np.count_nonzero(region_idx==int(atlas_segment_id))
        volume_dict_names[segment_name] = np.sum(eroded_atlas_vol == int(atlas_segment_id)) * np.power(0.02,3)
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

# Save region cell counts to region_cell_counts.csv
df = pd.DataFrame([corrected_count_dict,corrected_density_dict])
fname = 'region_cell_counts_filtered_' + str(cell_detection_filter['size'][0]) + 'px.csv'
fname = os.path.join('/jukebox/witten/Chris/data/clearmap2',fpath,fname)
df.to_csv(fname,index=False)
print('Saving cell detection results broken down by brain region to:\n/jukebox/witten/Chris/data/clearmap2/' + fpath + '/' + fname)
print()
