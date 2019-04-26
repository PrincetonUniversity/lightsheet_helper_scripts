#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:48:54 2019

@author: wanglab
"""

import os
from ClearMap.cluster.utils import load_kwargs
from ClearMap.cluster.par_tools import output_analysis

fld = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/processed"

brains = [os.path.join(fld, xx) for xx in os.listdir(fld)]

for brain in brains:
    params = load_kwargs(brain)    
    params["annotationfile"] = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_75um_erosion.tif"
    for k,v in params.iteritems():
        if isinstance(v, basestring):
            print v
            if v[:13] == "/home/wanglab": 
                v = "/jukebox" + v[13:]
                params[k] = v
            elif v[:38] == "/mnt/bucket/PNI-centers/LightSheetData":
                v = "/jukebox" + v[23:]
                params[k] = v
                
    output_analysis(threshold = (20, 900), row = (3,3), check_cell_detection = False, **params)