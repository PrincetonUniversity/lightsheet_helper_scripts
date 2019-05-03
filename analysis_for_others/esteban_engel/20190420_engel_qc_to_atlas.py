#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 13:55:49 2019

@author: wanglab
"""

import os, numpy as np, time, cv2
from skimage.external import tifffile
from collections import Counter
from scipy.ndimage import grey_dilation
#%%
##########################################################RUNS IN PYTHON 3##########################################################################
if __name__ == "__main__":

    #set up
    dst = "/jukebox/LightSheetData/pni_viral_vector_core/promoter_exp_qc"
    if not os.path.exists(dst): os.mkdir(dst)
    #make logs dir
    if not os.path.exists((os.path.join(dst, "errors"))): os.mkdir((os.path.join(dst, "errors")))
    src = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/processed"
    lst = [os.path.join(src, xx) for xx in os.listdir(src)]
    transform = "all";#both for regwatlas, and only affine for sig adn reg #"all", "single": don"t consider reg with sig at all
    verbose = True
    atl_pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/average_template_25_sagittal_forDVscans_z_thru_240.tif"
    alpha = 0.6
    #loop
    for fld in lst:
        start = time.time()
        print(fld)
        
        #####check cell transform
        #clearmap transformed to atlas cell dataframe
        try:
            dataframe = np.load(os.path.join(fld, "clearmap_cluster_output/cells_transformed_to_Atlas.npy"))
        
            #load and convert to single voxel loc
            zyx = np.asarray([str((int(xx[2]), int(xx[1]), int(xx[0]))) for xx in dataframe])
            zyx_cnt = Counter(zyx)
                                
            #map to atlas
            atl = tifffile.imread(atl_pth)
            atl_cnn = np.zeros_like(atl)
            errors = []
            for zyx,v in zyx_cnt.items():
                z,y,x = [int(xx) for xx in zyx.replace("(","",).replace(")","").split(",")]
                try:
                    atl_cnn[z,y,x] = v*100
                except Exception as e:
                    print(e)
                    errors.append(e)
            if len(errors)>0:
                with open(os.path.join(dst, "errors/{}_errors.txt".format(os.path.basename(fld))), "a") as flll:
                    for err in errors:
                        flll.write(str(err)+"\n")
                    flll.close()
            merged = np.stack([atl, atl_cnn, np.zeros_like(atl)], -1)
            
            #reorient to horizontal
            merged = (np.swapaxes(merged, 0, 2)) #alternatively could save this as an RBG image
            
            #opencv
            output = merged[...,0].copy() #red is atlas (axis 1)
            for z in range(merged[...,0].shape[0]):
                cells = grey_dilation(merged[z, :, :, 1], size=(3, 3))
                cv2.addWeighted(cells, alpha, output[z, :, :], 1-alpha, 0, output[z, :, :]) #green is cells (axis 1)
                
            #makes grayscale image of cells overlaid on atlas
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_atlas.tif".format(os.path.basename(fld))), output)      
            
            print("\ntook {} seconds to make merged maps\n".format(time.time()-start))
            
        except:
            print("\n tranformed cells do not exist, check if all jobs have run \n")