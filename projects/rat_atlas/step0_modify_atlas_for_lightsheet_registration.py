#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:54:48 2019

@author: wanglab
"""

import SimpleITK as sitk, numpy as np, tifffile, matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":

    #range to clip y
    yrng = (76, 822)
    xrng = (20, 475)
    zrng = (62, 364)

    #load in wh colored atlas and reorient
    orgn = "/jukebox/LightSheetData/brodyatlas/atlas/original/WHS_SD_rat_atlas_v3.nii"
    ann = np.fliplr(sitk.GetArrayFromImage(sitk.ReadImage(orgn)))
    #here we are reorienting to A-P, D-V orientation (how our images are taken) and THEN CROP
    ann = ann[::-1][zrng[0]:zrng[1],yrng[0]:yrng[1],xrng[0]:xrng[1]]
    #set destination
    nann = "/jukebox/LightSheetData/brodyatlas/atlas/modified/WHS_SD_rat_atlas_v3_anterior_up_DV.tif"
    tifffile.imsave(nann, ann)

    #load in composite brain
    org = "/jukebox/LightSheetData/brodyatlas/atlas/original/WHS_SD_rat_T2star_v1.01.nii"
    #apply as above
    auto = np.fliplr(sitk.GetArrayFromImage(sitk.ReadImage(org)))
    auto = auto[::-1][zrng[0]:zrng[1],yrng[0]:yrng[1],xrng[0]:xrng[1]]
    #set dsetination
    nauto = "/jukebox/LightSheetData/brodyatlas/atlas/modified/WHS_SD_rat_T2star_v1.01_anterior_up_DV.tif"
    tifffile.imsave(nauto, auto)

    #remove skull in atlas using mask
    from skimage.external import tifffile
    nann = "/jukebox/LightSheetData/brodyatlas/atlas/modified/WHS_SD_rat_atlas_v3_anterior_up_DV.tif"
    nauto = "/jukebox/LightSheetData/brodyatlas/atlas/modified/WHS_SD_rat_T2star_v1.01_anterior_up_DV.tif"
    mask = tifffile.imread(nann)
    atlas = tifffile.imread(nauto)
    atlas[mask==0] = 0
    skrm = "/jukebox/LightSheetData/brodyatlas/atlas/modified/WHS_SD_rat_T2star_v1.01_anterior_up_skullremoved_DV.tif"
    tifffile.imsave(skrm, atlas)

    #make df for labels
    ndf = pd.DataFrame(columns=["name", "id"])

    #wh labels - index = pix value
    wh ="/jukebox/LightSheetData/brodyatlas/atlas/original/WHS_SD_rat_atlas_v3.label"
    with open(wh, "r") as w:
        lines = w.readlines()
        w.close()

    #generate dataframe
    lines=lines[14:]
    for i in range(len(lines)):
        print ("{} of {}".format(i, len(lines)))
        line = lines[i]
        vals, structure, empty = line.split("'")
        idx, r, g, b, a, vis, msh = vals.split(); r=int(r); g=int(g); b=int(b)
        ndf.loc[i] = [structure, idx]

    ndf.to_csv("/jukebox/LightSheetData/brodyatlas/atlas/modified/labels_v3.csv")
