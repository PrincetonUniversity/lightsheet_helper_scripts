# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 23:02:31 2020

@author: Zahra
"""

#NOTE THIS ESSENTIALLY SCALES PIXEL SPACE*****
import numpy as np, os
from skimage.external import tifffile

#parallelized
jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

ann_path = '/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans_16bit.tif'
new_erode_path = '/jukebox/wang/zahra/kelly_cell_detection_analysis/eroded_atlases'
#
#ann_path = r'Y:\atlas\allen_atlas\annotation_2017_25um_sagittal_forDVscans_16bit.tif'
#new_erode_path = r'V:\kelly_cell_detection_analysis\eroded_atlases'

ann = tifffile.imread(ann_path)

struct_vals = np.unique(ann)[1:]
struct_microns_to_erode = 80
zyx_scale = (25,25,25)

iid = struct_vals[jobid]
print(iid)
sann = np.copy(ann)
sann[sann!=iid] = 0

from scipy.ndimage.morphology import distance_transform_edt
distance_space_inside = distance_transform_edt(sann.astype('bool'), sampling=zyx_scale)*-1 #INSIDE
distance_space_inside = np.abs(distance_space_inside)
mask = np.copy(distance_space_inside)
mask[distance_space_inside<=struct_microns_to_erode] = 0

#zero out edges
eann = np.copy(ann)
eann[mask==0]=0

tifffile.imsave(os.path.join(new_erode_path, "annotation_allen_2017_25um_sagittal_erode_80um_id{}.tif".format(iid)), eann.astype("uint16"))

print("\nannotation file saved for structure id : {}!".format(iid))