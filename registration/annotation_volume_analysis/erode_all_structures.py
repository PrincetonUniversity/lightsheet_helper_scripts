# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 23:02:31 2020

@author: Zahra
"""

#NOTE THIS ESSENTIALLY SCALES PIXEL SPACE*****
# %matplotlib inline
import numpy as np, os, matplotlib.pyplot as plt
from skimage.external import tifffile

#parallelized
jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

ann_path = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
new_erode_path = "/jukebox/wang/zahra/registration_error_pma/eroded_atlases"

#make fld
if not os.path.exists(new_erode_path): os.mkdir(new_erode_path)

# ann_path = r'Y:\atlas\allen_atlas\annotation_2017_25um_sagittal_forDVscans_16bit.tif'
# new_erode_path = r'Z:\zahra\kelly_cell_detection_analysis\eroded_atlases'

ann = tifffile.imread(ann_path)

struct_vals = np.unique(ann)
struct_microns_to_erode = 80
zyx_scale = (20,20,20)

iid = struct_vals[jobid]
print(iid)
sann = np.copy(ann)
sann[sann!=iid] = 0

#check
# plt.imshow(sann[300])

from scipy.ndimage.morphology import distance_transform_edt
distance_space_inside = distance_transform_edt(sann.astype('bool'), sampling=zyx_scale) #INSIDE

#check
# plt.imshow(distance_space_inside[300])

mask = np.copy(distance_space_inside)
mask[distance_space_inside<=struct_microns_to_erode] = 0


tifffile.imsave(os.path.join(new_erode_path, "annotation_pma_2018_20um_sagittal_erode_80um_id%d.tif" % iid), 
                mask.astype("bool").astype("float32")*iid)

print("\nannotation file saved for structure id : {}!".format(iid))