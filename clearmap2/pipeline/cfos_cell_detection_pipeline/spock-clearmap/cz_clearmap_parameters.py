# General imports
import pickle,sys

# ClearMap2 imports
sys.path.append('/jukebox/witten/Chris/python/ClearMap2-master')
import ClearMap.ImageProcessing.Experts.Cells as cells
import ClearMap.Utils.HierarchicalDict as hdict

# Set ClearMap2 cell detection parameters
cell_detection_parameter = cells.default_cell_detection_parameter.copy()
cell_detection_parameter['iullumination_correction'] = None
cell_detection_parameter['background_correction']['shape'] = (20,20)
cell_detection_parameter['background_correction']['form'] = 'Disk'
cell_detection_parameter['background_correction']['save'] = False
cell_detection_parameter['equalization'] = None
cell_detection_parameter['dog_filter'] = None
cell_detection_parameter['maxima_detection']['h_max'] = None
cell_detection_parameter['maxima_detection']['shape'] = 10
cell_detection_parameter['maxima_detection']['threshold'] = None
cell_detection_parameter['maxima_detection']['valid'] = True
cell_detection_parameter['maxima_detection']['save'] = False
cell_detection_parameter['shape_detection']['threshold'] = 750
cell_detection_parameter['shape_detection']['save'] = False
cell_detection_parameter['intensity_detection']['method'] = 'max'
cell_detection_parameter['intensity_detection']['shape'] = 10
cell_detection_parameter['intensity_detection']['measure'] = ['source','background']
cell_detection_parameter['verbose'] = False

# Save cell_detection_parameter.p
fname = '/jukebox/witten/Chris/data/clearmap2/utilities/cell_detection_parameter.p'
with open(fname, 'wb') as f:
    pickle.dump(cell_detection_parameter,f)
cell_detection_parameter = None

# Load and print cell_detection_parameter.p
with open(fname,'rb') as f:
    cell_detection_parameter = pickle.load(f)
print()
print('path                    : ' + fname)
hdict.pprint(cell_detection_parameter)
print()

# Set ClearMap2 cell detection filter thresholds
cell_detection_filter = {}
cell_detection_filter['x'] = (None,None)
cell_detection_filter['y'] = (None,None)
cell_detection_filter['z'] = (None,None)
cell_detection_filter['size'] = (20,None)
cell_detection_filter['intensity'] = (None,None)
cell_detection_filter['region'] = (None,None)

# Save cell_detection_filters.p
fname = '/jukebox/witten/Chris/data/clearmap2/utilities/cell_detection_filter.p'
with open(fname, 'wb') as f:
    pickle.dump(cell_detection_filter,f)
cell_detection_filter = None

# Load and print cell_detection_filters.p
with open(fname,'rb') as f:
    cell_detection_filter = pickle.load(f)
print()
print('path     : ' + fname)
hdict.pprint(cell_detection_filter)
print()
