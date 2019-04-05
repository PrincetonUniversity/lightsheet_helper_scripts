#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 17:35:02 2019

@author: wanglab
"""

import os, sys
import cPickle as pickle
os.chdir("/jukebox/wang/zahra/clearmap_cluster_copy")
from ClearMap.cluster.utils import load_kwargs
from ClearMap.ImageProcessing.CellDetection import detectCells
from ClearMap.cluster.preprocessing import pth_update
from ClearMap.parameter_file import set_parameters_for_clearmap
from ClearMap.cluster.par_tools import join_results_from_cluster, output_analysis

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
pth = "/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/processed"

flds = os.listdir(pth)

for fld in flds:
    brain = os.path.join(pth, fld)
    params = load_kwargs(brain)
    print(brain)
    
    #changing params from cluster - not ideal
    params["packagedirectory"] = os.path.join(brain, "clearmap_cluster")
    params["parameterfolder"] = os.path.join(brain, "clearmap_cluster/parameterfolder")
    
    #save params
    save_kwargs(**params)
    
    #changing cell detection param
    dct = pth_update(set_parameters_for_clearmap(testing = False, **params))
    dct["ImageProcessingParameter"]["detectSpotsParameter"]["removeBackgroundParameter"]["size"] = (5, 5)
    dct["ImageProcessingParameter"]["detectSpotsParameter"]["detectCellShapeParameter"]["threshold"] = 350
    
#########################################################################STEP 4##############################################################################
    
    for jobid in range(20): #assuming it doesnt use more than 20 chunks
        dct["ImageProcessingParameter"]["jobid"]=jobid
        #detect cells
        try: 
            result, substack = detectCells(**dct["ImageProcessingParameter"])
            if result == "ENDPROCESS": print("Jobid > # of jobs required, ending job")
        except:
            print("Jobid > # of jobs required, ending job")
    
    print("\n           finished step 4 - cell detection \n")                                     
#########################################################################STEP 5##############################################################################
    
    join_results_from_cluster(**params)    
    print("\n           finished step 5 \n")   
#########################################################################STEP 6##############################################################################

    try:
        output_analysis(threshold = (20, 900), row = (3,3), check_cell_detection = False, **params)
        print("\n           finished step 6 \n")   
    except:
        print("\n           using PMA atlas which is not currently compatible with clearmap functions. run in group separately")
                                