#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:48:54 2019

@author: wanglab
"""

import os, sys
sys.path.append("/jukebox/wang/zahra/clearmap_cluster_copy")
from ClearMap.cluster.utils import load_kwargs
from ClearMap.cluster.par_tools import output_analysis
import cPickle as pickle

def save_kwargs(pckloc = None, verbose = False, **kwargs):
    """
    Save out kwargs as param_dict.p unless otherwise noted.

    Inputs
    ----------------
    pckloc (optional) = location and name to save .p file
    kwargs = dictionary to pickle

    Returns
    ----------------
    location/name of file
    """
    #handle input
    if not pckloc: pckloc=os.path.join(kwargs["outputdirectory"], "param_dict.p")

    #verbosity
    if verbose: sys.stdout.write("\n\npckloc: {}".format(pckloc))

    #pickle
    pckfl=open(pckloc, "wb"); pickle.dump(kwargs, pckfl); pckfl.close()

    return pckloc

#%%

#NOTE: I COMMENTED OUT THE PART WHERE CELLS ARE TRANSFORMED IN CLEARMAP - DON'T NEED THAT    
fld = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/processed"

brains = [os.path.join(fld, xx) for xx in os.listdir(fld) if not xx == "v145_none_c3_1"]

for brain in brains:
    
    params = load_kwargs(brain)    
    params["annotationfile"] = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_75um_erosion_100um_ventricular_erosion.tif"
                
    save_kwargs(verbose = True, **params)
                
    output_analysis(threshold = (20, 900), row = (3,3), check_cell_detection = False, **params)
    
#%%    
brain = os.path.join(fld, "v145_none_c3_1")

params = load_kwargs(brain)    
params["annotationfile"] = "/jukebox/LightSheetData/pni_viral_vector_core/201808_promoter_exp/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_zflipped_75um_erosion_100um_ventricular_erosion.tif"
            
save_kwargs(verbose = True, **params)
            
output_analysis(threshold = (20, 900), row = (3,3), check_cell_detection = False, **params)

#%%

    
fld = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/processed"

brains = [os.path.join(fld, xx) for xx in os.listdir(fld)]

for brain in brains:
    
    params = load_kwargs(brain)    
    params["annotationfile"] = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans_z_thru_240_75um_erosion_100um_ventricular_erosion.tif"
                
    save_kwargs(verbose = True, **params)
                
    output_analysis(threshold = (20, 900), row = (3,3), check_cell_detection = False, **params)
    