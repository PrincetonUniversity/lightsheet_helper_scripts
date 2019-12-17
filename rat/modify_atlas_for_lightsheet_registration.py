#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:54:48 2019

@author: wanglab
"""

import nibabel as nib, matplotlib.pyplot as plt, tifffile, numpy as np, pandas as pd

if __name__ == "__main__":
    
    dst = '/jukebox/LightSheetData/brodyatlas/atlas/modified'
    ann_pth = '/jukebox/LightSheetData/brodyatlas/atlas/original/WHS_SD_rat_atlas_v3.nii'
    atl_pth = '/jukebox/LightSheetData/brodyatlas/atlas/original/WHS_SD_rat_T2star_v1.01.nii'
    
    ann = nib.load(ann_pth)
    
    ann = ann.get_fdata()
    
    plt.imshow(np.swapaxes(np.flip(np.rot90(ann, 2), 2), 0,2)[250]) #horizontal
    
    #makes anterior up orientation of horizontal brain
    tifffile.imsave(dst+'/WHS_SD_rat_atlas_v3.tif', np.swapaxes(np.flip(np.rot90(ann, 2), 2), 0,2).astype("uint16"))
    
    #range to clip x,y
    xrng = (20, 475)
    yrng = (76, 850) #ranges based on 2018 atlas optimization, see lightsheet_helper_scripts/registration/rat_registration_to_waxholm_mri.py
    zrng = (62, 364)
    
    ann = np.swapaxes(np.flip(np.rot90(ann, 2), 2), 0,2)
    ann = ann[zrng[0]:zrng[1],yrng[0]:yrng[1],xrng[0]:xrng[1]]
    tifffile.imsave(dst+'/WHS_SD_rat_atlas_v3_cropped.tif', ann.astype("uint16"))
    
    #convert to sagittal
    plt.imshow(np.swapaxes(ann, 0,2)[250]) #sagittal
    tifffile.imsave(dst+'/WHS_SD_rat_atlas_v3_cropped_sagittal.tif', np.swapaxes(ann, 0,2).astype("uint16"))
        
    #make df for labels
    ndf = pd.DataFrame(columns=["name", "id"])
    
    #wh labels - index = pix value
    wh ='/jukebox/LightSheetData/brodyatlas/atlas/original/WHS_SD_rat_atlas_v3.label'
    with open(wh, 'r') as w:
        lines = w.readlines()
        w.close()
    
    #generate dataframe
    lines=lines[14:]
    for i in range(len(lines)):
        print ('{} of {}'.format(i, len(lines)))
        line = lines[i]
        vals, structure, empty = line.split('"')
        idx, r, g, b, a, vis, msh = vals.split(); r=int(r); g=int(g); b=int(b)
        ndf.loc[i] = [structure, idx]
    
    ndf.to_csv('/jukebox/LightSheetData/brodyatlas/atlas/modified/labels_v3.csv')
    
