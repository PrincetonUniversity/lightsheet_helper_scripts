#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 17:35:02 2019

@author: wanglab
"""

import os, sys, pickle
from itertools import product
sys.path.append("/jukebox/wang/zahra/python/ClearMapCluster")
from ClearMap.cluster.utils import load_kwargs
from ClearMap.ImageProcessing.CellDetection import detectCells
from ClearMap.cluster.preprocessing import pth_update
from ClearMap.parameter_file import set_parameters_for_clearmap
from ClearMap.cluster.par_tools import join_results_from_cluster_helper, output_analysis
import ClearMap.IO as io
from ClearMap.Analysis.Statistics import thresholdPoints

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

if __name__ == "__main__":
    
    #get job id from SLURM
    print(sys.argv)
    print(os.environ["SLURM_ARRAY_TASK_ID"])
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"]) #to iterate through parameters
    stepid = int(sys.argv[1]) #to iterature through brains (790 dorsal + ventral)
    
    #brains to iterate through
    flds = ["/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos/171209_f37104_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_15-46-13",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos/171210_m37079_mouse2_20171014_790_015na_1hfsds_z5um_1000msec_17-25-29",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos/171215_m37081_observer_20171014_790_015na_1hfsds_z5um_1000msec_11-18-43",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810/181011_f37077_observer_20171011_790_017na_1hfds_z5um_1000msec_13-29-49",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810/181011_f37104_demonstrator_20171016_790_017na_1hfds_z5um_1000msec_11-55-12",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810/181015_f37080_mouse1_20171015_790_017na_1hfds_z5um_1000msec_08-22-40"]
    
    #params to iterate through - refer to Reiner et al or core facility notes for more info
    rBP_size_r = [3,5,7]
    fEMP_hmax_r = [None]
    fEMP_size_r = [0]
    fEMP_threshold_r = [None]
    fIP_method_r = ["Max"] 
    fIP_size_r = [5,7,10,12,15]
    dCSP_threshold_r = [50,75,100]
    
    # calculate number of iterations
    tick = 0
    for rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold in product(rBP_size_r,fEMP_hmax_r, fEMP_size_r,
        fEMP_threshold_r, fIP_method_r, fIP_size_r, dCSP_threshold_r):
        tick +=1
    sys.stdout.write("\n\nNumber of iterations is {}:".format(tick))
        
    #select brain based on first command line arg
    fld = flds[stepid]
    brain = os.path.join(fld, "clearmap")
    #load kwargs
    params = load_kwargs(brain)
    print("\n************"+os.path.dirname(brain)+"************\n")
    #changing params from cluster - not ideal
    params["packagedirectory"] = os.path.join(brain, "ClearMapCluster")
    params["parameterfolder"] = os.path.join(brain, "ClearMapCluster/parameterfolder")
    #save params
    save_kwargs(**params)
    #make cell detection dictionary
    dct = pth_update(set_parameters_for_clearmap(testing = False, **params))
    
    #select params based on jobid
    rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold=[(rBP_size,
        fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold) for rBP_size, fEMP_hmax,
        fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold in product(rBP_size_r, fEMP_hmax_r, fEMP_size_r,
        fEMP_threshold_r, fIP_method_r, fIP_size_r, dCSP_threshold_r)][jobid]
    #changing cell detection params
    #set params for sweep
    kwargs["removeBackgroundParameter_size"] = (rBP_size,rBP_size) #Remove the background with morphological opening (optimised for spherical objects), e.g. (7,7)
    kwargs["findExtendedMaximaParameter_hmax"] = fEMP_hmax # (float or None)     h parameter (for instance 20) for the initial h-Max transform, if None, do not perform a h-max transform
    kwargs["findExtendedMaximaParameter_size"] = fEMP_size # size in pixels (x,y) for the structure element of the morphological opening
    kwargs["findExtendedMaximaParameter_threshold"] = fEMP_threshold # (float or None)     include only maxima larger than a threshold, if None keep all local maxima
    kwargs["findIntensityParameter_method"] =  fIP_method # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
    kwargs["findIntensityParameter_size"] = (fIP_size,fIP_size,fIP_size) # (tuple)             size of the search box on which to perform the *method*
    kwargs["detectCellShapeParameter_threshold"] = dCSP_threshold # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated
    
    dct["ImageProcessingParameter"]["detectSpotsParameter"]["detectCellShapeParameter"]["threshold"] = 150
    
#########################################STEP 4###########################################################
    for iterid in range(50): #assuming it doesnt use more than 20 chunks
        dct["ImageProcessingParameter"]["jobid"]=iterid
        #detect cells
        try: 
            result, substack = detectCells(**dct["ImageProcessingParameter"])
            if result == "ENDPROCESS": print("Jobid > # of jobs required, ending job")
        except:
            print("Jobid > # of jobs required, ending job")
    
    print("\n           finished step 4 - cell detection \n")                                     
###########################################STEP 5########################################################
    
    out = join_results_from_cluster_helper(**dct["ImageProcessingParameter"])    
    print("\n           finished step 5 \n")   
###########################################STEP 6########################################################
    
    threshold = (500,3000); row = (2,2)
    points, intensities = io.readPoints(dct["ImageProcessingParameter"]["sink"])
    #Thresholding: the threshold parameter is either intensity or size in voxel, depending on the chosen "row"
    #row = (0,0) : peak intensity from the raw data
    #row = (1,1) : peak intensity from the DoG filtered data
    #row = (2,2) : peak intensity from the background subtracted data
    #row = (3,3) : voxel size from the watershed
    points, intensities = thresholdPoints(points, intensities, threshold = threshold, row = row);
    #points, intensities = thresholdPoints(points, intensities, threshold = (20, 900), row = (2,2));
    io.writePoints(dct["FilteredCellsFile"], (points, intensities));
    print("\n           finished step 6 \n")   
                                    