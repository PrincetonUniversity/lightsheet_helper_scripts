#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 11:50:58 2019

@author: wanglab
"""

#Modifications of Tom"s clearmap cluster code to parameter sweep an processed brain locally

import os, sys, tifffile, numpy as np
from skimage.exposure import rescale_intensity
from itertools import product
os.chdir("/jukebox/wang/zahra/clearmap_cluster_copy")
from ClearMap.cluster.par_tools import celldetection_operations
from ClearMap.cluster.utils import load_kwargs
from ClearMap.cluster.preprocessing import makedir, listdirfull, removedir


#sweep parameters copy & modifications 
def sweep_parameters_cluster(dst, jobid, optimization_chunk=12, pth=False, rescale=False, cleanup=True, **kwargs):
    """Function to sweep parameters
    
    final outputs will be saved in outputdirectory/parameter_sweep
    second copy will be saved in outputdirectory/parameter_sweep_jobid if cleanup=False

    Inputs:
        ----------------
        jobid: chunk of tissue to run (usually int between 20-30)
        #pth (optional): if pth to output folder after running package, function will load the param file automatically
        rescale (optional): str of dtype to rescale to. E.g.: "uint8"
        cleanup = T/F removes subfolders after
        optimization_chunk = this was the old "jobid" in this case it is the chunk of volume to look at
        kwargs (if not pth): "params" from run_clearmap_cluster.py
    """

    ######################################################################################################
    #NOTE: To adjust parameter sweep, modify ranges below
    ######################################################################################################
#    #first - with cleanup=True
#    rBP_size_r = range(3,9,2) #[5, 11] #range(5,19,2) ###evens seem to not be good <-- IMPORTANT TO SWEEP
    fEMP_hmax_r = [None]#[None, 5, 10, 20, 40]
    fEMP_size_r = [5]#range(3,8)
    fEMP_threshold_r = [None] #range(0,10)
    fIP_method_r = ["Max"] #["Max, "Mean"]
    fIP_size_r = [5]#range(1,5)
#    dCSP_threshold_r = range(10, 105, 20) #<-- IMPORTANT TO SWEEP
#    
    #second cleanup=False
    rBP_size_r = [5] #zmd commented out
    dCSP_threshold_r = [35]
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    
    
    # calculate number of iterations
    tick = 0
    for rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold in product(rBP_size_r, fEMP_hmax_r, fEMP_size_r, fEMP_threshold_r, fIP_method_r, fIP_size_r, dCSP_threshold_r):
        tick +=1

    sys.stdout.write("\n\nNumber of iterations is {}:".format(tick))
    
    #if pth is set - zmd added
    if pth:
        kwargs = load_kwargs(pth)
        
    #make folder for final output:
    opt = dst; makedir(opt)
    out = opt+"/parameter_sweep"; makedir(out)
    out0 = opt+"/parameter_sweep_jobid_{}".format(str(jobid).zfill(4)); makedir(out0)

    ntick = 0
    rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold=[(rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, 
                    fIP_method, fIP_size, dCSP_threshold) for rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, 
                    dCSP_threshold in product(rBP_size_r, fEMP_hmax_r, fEMP_size_r, fEMP_threshold_r, fIP_method_r, fIP_size_r, dCSP_threshold_r)][jobid]

    #zmd modified
    pth = out0+"/parametersweep_rBP_size{}_fEMP_hmax{}_fEMP_size{}_fEMP_threshold{}_fIP_method{}_fIP_size{}_dCSP_threshold{}.tif".format(rBP_size, 
                                        fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold)

    if not os.path.exists(pth):

        try:

            #set params for sweep
            kwargs["removeBackgroundParameter_size"] = (rBP_size, rBP_size) #Remove the background with morphological opening (optimised for spherical objects), e.g. (7,7)
            kwargs["findExtendedMaximaParameter_hmax"] = fEMP_hmax # (float or None)     h parameter (for instance 20) for the initial h-Max transform, if None, do not perform a h-max transform
            kwargs["findExtendedMaximaParameter_size"] = fEMP_size # size in pixels (x,y) for the structure element of the morphological opening
            kwargs["findExtendedMaximaParameter_threshold"] = fEMP_threshold # (float or None)     include only maxima larger than a threshold, if None keep all local maxima
            kwargs["findIntensityParameter_method"] =  fIP_method # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
            kwargs["findIntensityParameter_size"] = (fIP_size,fIP_size,fIP_size) # (tuple)             size of the search box on which to perform the *method*
            kwargs["detectCellShapeParameter_threshold"] = dCSP_threshold # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated

            #tmp
            import cPickle as pickle
            nkwargs = load_kwargs(kwargs["outputdirectory"])
            kwargs["outputdirectory"] = out0
            nkwargs.update(kwargs)
            pckloc=out0+"/param_dict.p"; pckfl=open(pckloc, "wb"); pickle.dump(nkwargs, pckfl); pckfl.close()

            #run cell detection
            ntick+=1
            sys.stdout.write("\n\n\n           *****Iteration {} of {}*****\n\n\n".format(ntick, tick))
            sys.stdout.write("    Iteration parameters: {}     {}     {}     {}     {}     {}     {}".format(kwargs["removeBackgroundParameter_size"], 
                             kwargs["findExtendedMaximaParameter_hmax"], kwargs["findExtendedMaximaParameter_size"], kwargs["findExtendedMaximaParameter_threshold"],         kwargs["findIntensityParameter_method"],         kwargs["findIntensityParameter_size"],        kwargs["detectCellShapeParameter_threshold"]))
            celldetection_operations(optimization_chunk, testing = True, **kwargs)

            #list, load, and maxip
            if ntick == 1: raw = [xx for xx in listdirfull(out0+"/optimization/raw") if "~" not in xx and ".db" not in xx]; raw.sort()
            raw_im = np.squeeze(tifffile.imread(raw)); raw_mx = np.max(raw_im, axis = 0)
            bkg = [xx for xx in listdirfull(out0+"/optimization/background") if "~" not in xx and "Thumbs.db" not in xx]; bkg.sort()
            bkg_im = tifffile.imread(bkg); bkg_mx = np.max(bkg_im, axis = 0)
            cell = [xx for xx in listdirfull(out0+"/optimization/cell") if "~" not in xx and ".db" not in xx]; cell.sort()
            cell_im = tifffile.imread(cell); cell_mx = np.max(cell_im, axis = 0)

            #optional rescale:
            if rescale:
                raw_mx = rescale_intensity(raw_mx, in_range=str(raw_mx.dtype), out_range=rescale).astype(rescale)
                bkg_mx = rescale_intensity(bkg_mx, in_range=str(bkg_mx.dtype), out_range=rescale).astype(rescale)
                cell_mx = rescale_intensity(cell_mx, in_range=str(cell_mx.dtype), out_range=rescale).astype(rescale)


            #concatenate and save out:
            bigim = np.concatenate((raw_mx, bkg_mx, cell_mx), axis = 1); del bkg, bkg_im, bkg_mx, cell, cell_im,cell_mx
            if cleanup: removedir(out0)
            if not cleanup: tifffile.imsave(pth, bigim, compress = 1)
            
            #save in main
            npth = out+"/jobid_{}_parametersweep_rBP_size{}_fEMP_hmax{}_fEMP_size{}_fEMP_threshold{}_fIP_method{}_fIP_size{}_dCSP_threshold{}.tif".format(str(jobid).zfill(4), 
                               rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold)
            tifffile.imsave(npth, bigim, compress = 1)
            

        except Exception, e:
            print ("Error on: {}\n\nerror={}".format(pth,e))
            im = np.zeros((10,10,10))
            tifffile.imsave(pth, im, compress = 1)
            with open(os.path.join(out, "errored_files.txt"), "a") as fl:
                fl.write("\n\n{}\n{}\n".format(pth, kwargs))
                fl.close

    return

#%%

brain = "/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos/201701_mk06"
kwargs = load_kwargs(brain)
dst = "/jukebox/wang/pisano/tracing_output/cfos/201902_reim_201701_cfos_analysis/parameter_sweep/201701_mk06"
#parameter sweep cell detection parameters. NOTE read all of functions description before using. VERY CPU intensive
#for first pass at cell detection
for jobid in range(1): #to find the range value, run lines 152 - 173 in this script, part of sweep_parameters_cluster func
    
    sweep_parameters_cluster(dst, jobid, cleanup = False, **kwargs)
