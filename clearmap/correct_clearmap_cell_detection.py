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
    print(os.environ["SLURM_ARRAY_TASK_ID"])
    stepid = int(os.environ["SLURM_ARRAY_TASK_ID"]) #to iterate through brains (790 dorsal + ventral)
    
    #brains to iterate through
    flds = ["/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos/171209_f37104_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_15-46-13",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos/171210_m37079_mouse2_20171014_790_015na_1hfsds_z5um_1000msec_17-25-29",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos/171215_m37081_observer_20171014_790_015na_1hfsds_z5um_1000msec_11-18-43",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810/181011_f37077_observer_20171011_790_017na_1hfds_z5um_1000msec_13-29-49",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810/181011_f37104_demonstrator_20171016_790_017na_1hfds_z5um_1000msec_11-55-12",
            "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810/181015_f37080_mouse1_20171015_790_017na_1hfds_z5um_1000msec_08-22-40"]
    
    #params to iterate through - refer to Reiner et al or core facility notes for more info
    rBP_size_r = [3,5,7]
    fIP_method_r = ["Max"] 
    fIP_size_r = [5,10,15]
    dCSP_threshold_r = [50,75,100]
    
    # calculate number of iterations
    tick = 0
    for rBP_size, fIP_method, fIP_size, dCSP_threshold in product(rBP_size_r,
        fIP_method_r, fIP_size_r, dCSP_threshold_r):
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
    
    for jobid in range(tick): #iterate through each parameter set
        #select params based on jobid
        rBP_size, fIP_method, fIP_size, dCSP_threshold=[(rBP_size,
            fIP_method, fIP_size, dCSP_threshold) for rBP_size, 
            fIP_method, fIP_size, dCSP_threshold in product(rBP_size_r,
            fIP_method_r, fIP_size_r, dCSP_threshold_r)][jobid]
        
        #changing cell detection params
        dct["ImageProcessingParameter"]["detectSpotsParameter"]["removeBackgroundParameter"]["size"] = (rBP_size,rBP_size)
        dct["ImageProcessingParameter"]["detectSpotsParameter"]["findIntensityParameter"]["method"] = fIP_method
        dct["ImageProcessingParameter"]["detectSpotsParameter"]["findIntensityParameter"]["size"] = (fIP_size, fIP_size, fIP_size)
        dct["ImageProcessingParameter"]["detectSpotsParameter"]["detectCellShapeParameter"]["threshold"] = dCSP_threshold
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
        points, intensities = thresholdPoints(points, intensities, threshold = threshold, 
                            row = row)
        #change dst to match parameters sweeped
        dst = (os.path.join(brain, 
        "clearmap_cluster_output/cells_rBP%s_fIPmethod%s_fIPsize%s_dCSP%s.npy" % (rBP_size,
            fIP_method, fIP_size, dCSP_threshold)),
        os.path.join(brain, 
        "clearmap_cluster_output/intensities_rBP%s_fIPmethod%s_fIPsize%s_dCSP%s.npy" % (rBP_size,
            fIP_method, fIP_size, dCSP_threshold)))
        
        io.writePoints(dst, (points, intensities));
        print("\n           finished step 6 \n")   
                                    