#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 11:14:42 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl
from tools.registration.register import transformed_pnts_to_allen_helper_func, count_structure_lister
from tools.registration.register import change_transform_parameter_initial_transform
from tools.registration.transform_list_of_points import create_text_file_for_elastix, modify_transform_files
from tools.registration.transform_list_of_points import point_transformix, unpack_pnts
from tools.utils.io import makedir, load_kwargs
from skimage.external import tifffile
from tools.analysis.network_analysis import make_structure_objects
from scipy.ndimage.measurements import center_of_mass
import matplotlib.colors
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop
 
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

def find_site(im, thresh=3, filter_kernel=(3,3,3), num_sites_to_keep=1):
    """Find a connected area of high intensity, using a basic filter + threshold + connected components approach
    
    by: bdeverett
    Parameters
    ----------
    img : np.ndarray
        3D stack in which to find site (technically need not be 3D, so long as filter parameter is adjusted accordingly)
    thresh: float
        threshold for site-of-interest intensity, in number of standard deviations above the mean
    filter_kernel: tuple
        kernel for filtering of image before thresholding
    num_sites_to_keep: int, number of injection sites to keep, useful if multiple distinct sites
    
    Returns
    --------
    bool array of volume where coordinates where detected
    """
    from scipy.ndimage.filters import gaussian_filter as gfilt
    from scipy.ndimage import label
    if type(im) == str: im = tifffile.imread(im)

    filtered = gfilt(im, filter_kernel)
    thresholded = filtered > filtered.mean() + thresh*filtered.std() 
    labelled,nlab = label(thresholded)

    if nlab == 0:
        raise Exception('Site not detected, try a lower threshold?')
    elif nlab == 1:
        return labelled.astype(bool)
    elif num_sites_to_keep == 1:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        return labelled == np.argmax(sizes)+1
    else:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        vals = [i+1 for i in np.argsort(sizes)[-num_sites_to_keep:][::-1]]
        return np.in1d(labelled, vals).reshape(labelled.shape)
    
dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/rtn"

ann_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/atlases/sagittal_allen_ann_25um_iso_60um_edge_160um_ventricular_erosion.tif"

#cut annotation file in middle
ann = tifffile.imread(ann_pth)
plt.imshow(ann[300])
z,y,x = ann.shape
#make sure each halves are same dimension as original ann
ann_left = np.zeros_like(ann)
ann_left[:int(z/2), :, :] = ann[:int(z/2), :, :] #cut in the middle in x
ann_right = np.zeros_like(ann)
ann_right[int(z/2):, :, :] = ann[int(z/2):, :, :]
plt.imshow(ann_left[120])

#collect 
#brains should be in this order as they were saved in this order for inj analysis
brains = ["20160622_db_bl6_crii_52hr_01", "20160622_db_bl6_unk_01", "20160801_db_cri_02_1200rlow_52hr", 
          "20160801_db_l7_cri_01_mid_64hr", "20160822_tp_bl6_crii_250r_01", "20160822_tp_bl6_crii_1250r_05", 
          "20160822_tp_bl6_crii_1500r_06", "20160823_tp_bl6_cri_250r_01",
          "20160823_tp_bl6_cri_500r_02", "20160912_tp_bl6_lob7_750r_04", "20160916_tp_lob6_ml_01", "20160916_tp_lob6_ml_04", 
          "20160916_tp_lob7_250r_05", "20160920_tp_bl6_lob7_250r_02", "20160920_tp_bl6_lob7_500r_03", "20160920_tp_bl6_lob7_ml_01", 
          "20161201_db_bl6_lob6b_500r_53d5hr", "20161203_tp_bl6_crii_250r_08", 
          "20161203_tp_bl6_lob7_500r_05", "20161203_tp_bl6_lob7_1000r_06", "20161203_tp_bl6_lob7_1500r_07",
          "20161203_tp_bl6_lob7_ml_04", "20161205_tp_bl6_sim_250r_02", "20161205_tp_bl6_sim_750r_03",
          "20161205_tp_bl6_sim_1250r_04", "20161205_tp_bl6_sim_1750r_05","20161207_db_bl6_lob6a_500r_53hr", 
          "20161207_db_bl6_lob6a_850r_53hr", "20161208_db_bl6_cri_50r_pv_53hr","20170115_tp_bl6_lob6b_500r_05",
          "20170116_tp_bl6_lob6b_lpv_07","20170130_tp_bl6_sim_rlat_05",
          "20170207_db_bl6_crii_1300r_02", "20170212_tp_bl6_lob8_ml_05",
          "20170308_tp_bl6f_cri_2x_03","20170308_tp_bl6f_lob7_2x_02", 
          "20170410_tp_bl6_lob6a_ml_repro_01", "20170410_tp_bl6_lob6a_ml_repro_02",
          "20170411_db_bl6_crii_lat_53hr", "20170411_db_bl6_crii_mid_53hr",
          "20170411_db_bl6_crii_rpv_53hr", "20170419_db_lob6b_rpv_53hr"]    

src = "/jukebox/wang/pisano/tracing_output/antero_4x"

imgs = [os.path.join(src, xx) for xx in brains]

#pool brain names and L/R designation into dict
lr_dist = {}
thal_inj_vol = {}

#save inj segments in folder
sv_dst = os.path.join(dst, "injection"); makedir(sv_dst)
#get inj vol roundabout way
for img in imgs:
    brain = os.path.basename(img)
    print(brain)
    kwargs = load_kwargs(img)
    try:
        inj_vol_pth = [xx for xx in kwargs["volumes"] if xx.ch_type == "injch"][0]
        inj_vol = tifffile.imread(inj_vol_pth.ch_to_reg_to_atlas)
        inj_vol[inj_vol < 0] = 10000
        assert np.sum(inj_vol < 0) == 0
        z,y,x = inj_vol.shape
        
        arr = find_site(inj_vol[:, 423:, :])
        #save segment
        tifffile.imsave(os.path.join(sv_dst, brain+".tif"), arr.astype("uint16"))
        
        z_c,y_c,x_c = center_of_mass(arr)
        #take distance from center to arbitrary "midline" (aka half of z axis)
        dist = z_c-(z/2)
        #save to dict 
        lr_dist[brain[:-8]] = dist
        thal_inj_vol[brain[:-8]] = np.sum(inj_vol)
        
        if dist < 0:
            print("brain {} has a left-sided injection\n".format(brain))
        elif dist > 0:
            print("brain {} has a right-sided injection\n".format(brain))
        else:
            print("brain has an injection close to midline so not considering it rn\n")
    except:
        print("brain %s has no injection volume, segment from elsewhere\n" % brain)


#make structures
#FIXME: for some reason the allen table does not work on this, is it ok to use PMA...    
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"

structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
lr_brains = list(lr_dist.keys())
atl_dst = os.path.join(dst, "pma_to_aba"); makedir(atl_dst)
id_table = pd.read_excel(df_pth)
