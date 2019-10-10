#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:27:56 2019

@author: wanglab
"""

import tifffile, matplotlib.pyplot as plt, numpy as np, os, sys
import matplotlib.colors
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "Red"])
from skimage import filters
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3") #if you clone this repo yourself you can then set this path to your copy of the code fyi
from tools.imageprocessing.orientation import fix_orientation
from tools.conv_net.utils.io import read_roi_zip
from tools.registration.transform import transformed_pnts_to_allen_helper_func
from tools.registration.register import count_structure_lister
from tools.utils.io import load_kwargs

brain = "/jukebox/wang/oostland/lightsheet/m26"

brainname = os.path.basename(brain)
roipth = os.path.join(brain, "probe_track_points_coronal.zip")

rois = np.array([[int(yy) for yy in xx[0].replace(".roi", "").split("-")] for xx in
                            read_roi_zip(roipth, include_roi_name=True)])

pnts = np.array([[xx[0]-1, xx[1], xx[2]] for xx in rois]) #correct for fijis 1-based z numbering
pnts_sag = np.array([[xx[2], xx[0], xx[1]] for xx in pnts])

kwargs = load_kwargs(brain)
impth = [vol for vol in kwargs["volumes"] if vol.ch_type == "injch" or vol.ch_type == "cellch"][0].ch_to_reg_to_atlas

im = tifffile.imread(impth)
imcor = fix_orientation(im, ("2", "0", "1")) #points clicked in coronal vol

track = np.zeros_like(imcor)
size = (1,20,1) #size of dilatation in zyx
otsu_factor=1

for pnt in pnts:
    #print pnt
    vol = np.copy(imcor[np.max((pnt[0]-size[0],0)):pnt[0]+size[0], np.max((pnt[1]-size[1],0)):pnt[1]+size[1], 
                        np.max((pnt[2]-size[2],0)):pnt[2]+size[2]])*1.0
    v=filters.threshold_otsu(vol)/float(otsu_factor)
    vol[vol<v]=0
    vol[vol>=v]=1
    nvol = np.maximum(track[np.max((pnt[0]-size[0],0)):pnt[0]+size[0], np.max((pnt[1]-size[1],0)):pnt[1]+size[1], 
                            np.max((pnt[2]-size[2],0)):pnt[2]+size[2]], vol)
    track[np.max((pnt[0]-size[0],0)):pnt[0]+size[0], np.max((pnt[1]-size[1],0)):pnt[1]+size[1], 
          np.max((pnt[2]-size[2],0)):pnt[2]+size[2]]=nvol

merged = np.stack([imcor.astype("uint16"), track.astype("uint16"), np.zeros_like(imcor.astype("uint16"))], -1)

tifffile.imsave(os.path.join(brain, "%s_probe_track_overlay.tif" % brainname), merged.astype("uint16"))

#export coordinates
if os.path.exists(os.path.join(brain, "{}_allen_coordinates.txt".format(brainname))): os.remove(os.path.join(brain, "{}_allen_coordinates.txt".format(brainname)))
with open(os.path.join(brain, "{}_allen_coordinates.txt".format(brainname)), "a") as txt:
    txt.write("\nAtlas coordinates (zyx) in the saggital orientation:\n%s\n" % pnts_sag)
    txt.write("\nAtlas coordinates (zyx) in the coronal orientation:\n%s" % pnts)

            
#convert to structure
annotation_file = kwargs["annotationfile"]
ann = tifffile.imread(annotation_file)
zpnts, ypnts, xpnts = np.nonzero(fix_orientation(track, ("1", "2", "0"))) #make it back to sagittal for mapping
points = transformed_pnts_to_allen_helper_func([(zi,ypnts[i],xpnts[i]) for i, zi in enumerate(zpnts)], ann, order = "ZYX")    

#make dataframe
lut_path = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx" #corresponds with the atlas, change if changing atlas
df = count_structure_lister(lut_path, *points)
df.to_excel(os.path.join(brain, "%s_allen_structures.xlsx" % brainname))
