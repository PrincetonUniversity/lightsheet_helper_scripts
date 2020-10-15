#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 16:01:06 2020

@author: wanglab
"""

import os,tifffile as tif,numpy as np,matplotlib.pyplot as plt

#path to stitched and atlas folders
stitched="/home/wanglab/LightSheetTransfer/tp/PRV_50hr-019/Ex_642_Em_2/stitched/RES(7565x5731x3270)/108240/108240_112259"
ann="/home/wanglab/scratch/zmd/PRV_50hr-019/transformed_annotations/single_tifs"
dst="/home/wanglab/Desktop"

#planes to max project over?
zplns=[1901,1961]
ims=[os.path.join(stitched,xx) for xx in os.listdir(stitched)]
ims.sort()
#find ims within zplns
ims=ims[zplns[0]:zplns[1]]
#takes a while...
im=np.array([tif.imread(im) for im in ims])
#same for anns
anns=[os.path.join(ann,xx) for xx in os.listdir(ann)]
anns.sort()
anns=anns[zplns[0]:zplns[1]]
#takes a while...
an=np.array([tif.imread(a) for a in anns])
#save
tif.imsave(os.path.join(dst,"prv50hr_z%s-%s.tif" %(zplns[0],zplns[1])),np.max(im,axis=0).astype("uint16"),compress=6)
tif.imsave(os.path.join(dst,"prv50hr_ann_z%s-%s.tif" %(zplns[0],zplns[1])),np.max(an,axis=0).astype("float32"),compress=6)
