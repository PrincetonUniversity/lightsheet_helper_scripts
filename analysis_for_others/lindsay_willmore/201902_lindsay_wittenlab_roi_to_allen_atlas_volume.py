#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 13:14:39 2019

@author: wanglab
"""

import os, numpy as np, time
from skimage.external import tifffile
import matplotlib.pyplot as plt
import matplotlib.colors

os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from tools.conv_net.input.read_roi import read_roi_zip
from tools.registration.transform import transformed_pnts_to_allen_helper_func
from tools.registration.register import count_structure_lister

if __name__ == "__main__":
    
    src = "/jukebox/wang/willmore/lightsheet/20190304_fiber_placement"
    
    flds = [os.path.join(src, xx) for xx in os.listdir(src) if os.path.isdir(os.path.join(src, xx))]
    
    for fld in flds:
        start = time.time()
          
        brain = os.path.basename(fld)
        #load in ROIS - clicked in horizontal volume
        roi_pth = fld + "/{}_20190314_fiber_points_RoiSet.zip".format(brain)
        #check if path exists obv
        if os.path.exists(roi_pth):
                
            zyx_rois = np.asarray([[int(yy) for yy in xx.replace(".roi", "").split("-")] for xx in read_roi_zip(roi_pth, include_roi_name=True)])
                
            #make merged image
            zyx = np.asarray([str((int(xx[0]), int(xx[1]), int(xx[2]))) for xx in zyx_rois])
            
            #make destination path
            dst = os.path.join(fld, "points_merged_to_atlas")
            if not os.path.exists(dst): os.mkdir(dst)
            
            #export coordinates
            if os.path.exists(os.path.join(dst, "{}_allen_coordinates.txt".format(brain))): os.remove(os.path.join(dst, "{}_allen_coordinates.txt".format(brain)))
            with open(os.path.join(dst, "{}_allen_coordinates.txt".format(brain)), "a") as txt:
                txt.write("Allen Atlas CCF coordinates (in zyx):\n%s" % zyx_rois)
            
            #atlas (horizontal)
            atl = tifffile.imread("/jukebox/wang/zahra/atlas/average_template_25_sagittal_forDVscans.tif")
            
            atl_cnn = np.zeros_like(atl)
            
            for i in range(zyx_rois.shape[0]):
                atl_cnn[zyx_rois[i][0], zyx_rois[i][1], zyx_rois[i][2]] = 1
         
            merged = np.stack([atl, atl_cnn, np.zeros_like(atl)], -1)
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_horizontal.tif".format(brain)), merged)
            
            coronal = np.transpose(merged, [1, 0, 2, 3]) #make coronal sections
            
            #next, make smaller sections to visualise site better
            z = np.nonzero(coronal[..., 1])[0]
            
            #find z range of sites
            site1 = z[0:5]; site2 = z[5:]
            zrange_site1 = range(min(site1)-1, max(site1)+2); zrange_site2 = range(min(site2)-1, max(site2)+2)
            
            #save out coronal sections - based on the fact that you click 5 points in each site in HORIZONTAL sectionss
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_coronal.tif".format(brain)), coronal)
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_coronal_site1_z{}_{}.tif".format(brain, min(zrange_site1), max(zrange_site1))), coronal[zrange_site1])
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_coronal_site2_z{}_{}.tif".format(brain, min(zrange_site2), max(zrange_site2))), coronal[zrange_site2])
            
            #doing a max projection, in case you just want to look at that
            maxip1 = np.max(coronal[zrange_site1], 0)
            maxip2 = np.max(coronal[zrange_site2], 0)
            
            alpha = 0.6 #determines transparency, don't need to alter
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #color makes cells pop
            plt.imshow(maxip1[...,0], "gist_yarg")
            plt.imshow(maxip1[...,1], cmap, alpha = 0.6)
            plt.savefig(os.path.join(dst, "{}_points_merged_to_Allen_coronal_site1_maxip_z{}_{}.pdf".format(brain, min(zrange_site1), max(zrange_site1))), dpi = 300)
            
            plt.imshow(maxip2[...,0], "gist_yarg")
            plt.imshow(maxip2[...,1], cmap, alpha = 0.6)
            plt.savefig(os.path.join(dst, "{}_points_merged_to_Allen_coronal_site2_maxip_z{}_{}.pdf".format(brain, min(zrange_site2), max(zrange_site2))), dpi = 300)
            
            print("\n\ntook {} seconds to make merged maps for {}\n".format(time.time()-start, brain))
            
            #make allen structure LUT
            #go from horiztonal to sag
            zyx_rois = np.asarray([[xx[2], xx[1], xx[0]] for xx in zyx_rois])
            
            #convert to structure
            annotation_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif"
            ann = tifffile.imread(annotation_file)
            points = transformed_pnts_to_allen_helper_func(list(zyx_rois), ann, order = "ZYX")    
            
            #make dataframe
            lut_path = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table.xlsx"
            df = count_structure_lister(lut_path, *points)
            df.to_excel(os.path.join(dst, "{}_allen_structures.xlsx".format(brain)))