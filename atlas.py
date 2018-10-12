#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:09:03 2018

@author: wanglab
"""

import numpy as np
from skimage.external import tifffile

#load atlas and annotations file
allen_ann = tifffile.imread('/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif') #load annotations

allen_atlas = tifffile.imread('/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif') #load atlas

#find atlas and annotations file shape - triple of z,y,x coordinates
ann_shape = allen_ann.shape

atlas_shape = allen_atlas.shape

#straight cut
#'zero' out cerebellum that doesn't exist in the scans
#allen_ann[:,390:,:] = np.zeros((456, 138, 320))
#
#allen_atlas[:,390:,:] = np.zeros((456, 138, 320))

#tilted cut - 201810 cfos data
#messy but trying to think it through
y1 = 461; y2 = 411; x1 = 0; x2 = 115

##annotations
#zeroed_matrix = allen_ann[:, y2:y1, x1:x2] #z,y,x convention
#for i in range(ann_shape[0]):
#    zeroed_matrix_lower_triangle = np.flip(np.tril(np.flip(zeroed_matrix[i,:,:],axis=0)),axis=0)
#    allen_ann[i, y2:y1, x1:x2] = zeroed_matrix_lower_triangle
##    allen_ann[i, y1:ann_shape[1], :] = np.zeros((1, ann_shape[1]-y1, ann_shape[2]))

#atlas
zeroed_matrix = allen_atlas[:, y2:y1, x1:x2] #z,y,x convention
for i in range(ann_shape[0]):
    zeroed_matrix_lower_triangle = np.flip(np.tril(np.flip(zeroed_matrix[i,:,:], axis=0)),axis=0)
    allen_atlas[i, y2:y1, x1:x2] = zeroed_matrix_lower_triangle
#    allen_atlas[i, y1:atlas_shape[1], :] = np.zeros((1, atlas_shape[1]-y1, atlas_shape[2]))

#save out tifs in the appropriate location
#tifffile.imsave('/jukebox/LightSheetTransfer/atlas/201810_jess_cfos/annotation_template_25_sagittal_forDVscans_forebrain_only_tilted_cut.tif', allen_ann)

tifffile.imsave('/jukebox/LightSheetTransfer/atlas/201810_jess_cfos/average_template_25_sagittal_forDVscans_forebrain_only_tilted_cut.tif', allen_atlas)