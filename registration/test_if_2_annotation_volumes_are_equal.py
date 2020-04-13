#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:38:46 2020

@author: wanglab
"""

import numpy as np, SimpleITK as sitk

p1 = '/home/wanglab/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans_do_not_use.nrrd'
p2 = '/home/wanglab/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans.tif'

ann1 = sitk.GetArrayFromImage(sitk.ReadImage(p1))
ann2 = sitk.GetArrayFromImage(sitk.ReadImage(p2))

#test to see if both volumes are unequal at any point
test = ann1!=ann2
#if the sum of the test array is 0, the 2 volumes are equal
assert np.sum(test) == 0

iids1 = np.unique(ann1).astype(int)
iids2 = np.unique(ann2).astype(int)
#check if they have the same annotations
assert np.sum(iids1!=iids2) == 0
