#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 13:14:39 2019

@author: wanglab
"""

import os, numpy as np, time
from skimage.external import tifffile
from collections import Counter
os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from tools.conv_net.input.read_roi import read_roi_zip
from tools.registration.register import elastix_command_line_call, jacobian_command_line_call, change_interpolation_order, transformix_command_line_call, count_structure_lister
from tools.registration.transform import transformed_pnts_to_allen_helper_func

if __name__ == "__main__":
    
    start = time.time()
    
    fld = "/jukebox/wang/willmore/lightsheet/20181213_fiber_optic_placement/processed/m784"
    
    #load in ROIS - clicked in horizontal volume
    roi_pth = "/jukebox/wang/willmore/lightsheet/20181213_fiber_optic_placement/processed/m784/20190211_fiber_points_horizontal_RoiSet.zip"
    zyx_rois = np.asarray([[int(yy) for yy in xx.replace(".roi", "").split("-")] for xx in read_roi_zip(roi_pth, include_roi_name=True)])
        
    #make merged image
    zyx = np.asarray([str((int(xx[0]), int(xx[1]), int(xx[2]))) for xx in zyx_rois])
    zyx_cnt = Counter(zyx)
    #atlas
    atl = tifffile.imread("/jukebox/wang/zahra/atlas/average_template_25_sagittal_forDVscans.tif")
    
    atl_cnn = np.zeros_like(atl)
    errors = []
    
    for i in range(zyx_rois.shape[0]):
        atl_cnn[zyx_rois[i][0], zyx_rois[i][1], zyx_rois[i][2]] = 100
        
    merged = np.stack([atl, atl_cnn, np.zeros_like(atl)], -1)
    tifffile.imsave(os.path.join(fld, "points_merged_to_atlas/{}_points_merged_to_Allen_horizontal.tif".format(os.path.basename(fld))), merged)
    
    print("\ntook {} seconds to make merged maps\n".format(time.time()-start))