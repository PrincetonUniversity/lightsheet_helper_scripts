# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 12:27:38 2020

@author: Zahra
"""

#NOTE THIS ESSENTIALLY SCALES PIXEL SPACE*****
#%matplotlib inline
import numpy as np, os, matplotlib.pyplot as plt
from skimage.external import tifffile

#parallelized
jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

ann_path = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans_16bit.tif"
new_erode_path = "/jukebox/wang/zahra/kelly_cell_detection_analysis/dilated_atlases"
#
#ann_path = r'Y:\atlas\allen_atlas\annotation_2017_25um_sagittal_forDVscans_16bit.tif'
#new_erode_path = r'Z:\zahra\kelly_cell_detection_analysis\dilated_atlases'

ann = tifffile.imread(ann_path)

struct_vals = np.unique(ann)[1:]
struct_microns_to_dilate = 80
zyx_scale = (25,25,25)

#parallelized
iid = struct_vals[jobid]
    
print(iid)
sann = np.copy(ann)
sann[sann!=iid] = 0

from scipy.ndimage.morphology import distance_transform_edt
distance_space_outside = distance_transform_edt(np.logical_not(sann.astype("bool")), sampling=zyx_scale) #INSIDE

mask = np.copy(distance_space_outside)
mask[distance_space_outside >= struct_microns_to_dilate] = 0

#zero out
eann = np.copy(ann)
eann[mask == 0] = 0

sum_ann = sann.astype("bool")+eann.astype("bool")

tifffile.imsave(os.path.join(new_erode_path, "annotation_allen_2017_25um_sagittal_erode_80um_id%d.tif" % iid), sum_ann.astype("uint16")*iid)

print("\nannotation file saved for structure id : {}!".format(iid))