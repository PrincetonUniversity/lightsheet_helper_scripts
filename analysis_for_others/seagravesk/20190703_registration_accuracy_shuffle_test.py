#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 14:58:23 2019

@author: wanglab

registers 2 resampled for elastix volumes 100 times each
can be used to calculate registration error and error associated with each allen brain structure
"""

import os, sys, time, subprocess as sp

def elastix_command_line_call(fx, mv, out, parameters, fx_mask=False, verbose=False):
    """Wrapper Function to call elastix using the commandline, this can be time consuming
    
    Inputs
    -------------------
    fx = fixed path (usually Atlas for "normal" noninverse transforms)
    mv = moving path (usually volume to register for "normal" noninverse transforms)
    out = folder to save file
    parameters = list of paths to parameter files IN ORDER THEY SHOULD BE APPLIED
    fx_mask= (optional) mask path if desired
    
    Outputs
    --------------
    ElastixResultFile = ".tif" or ".mhd" result file
    TransformParameterFile = file storing transform parameters
    
    """
    e_params=["elastix", "-f", fx, "-m", mv, "-out", out]
    if fx_mask: e_params=["elastix", "-f", fx, "-m", mv, "-fMask", fx_mask, "-out", out]
    
    ###adding elastix parameter files to command line call
    for x in range(len(parameters)):
        e_params.append("-p")
        e_params.append(parameters[x])
    
    #set paths
    TransformParameterFile = os.path.join(out, "TransformParameters.{}.txt".format((len(parameters)-1)))
    ElastixResultFile = os.path.join(out, "result.{}.tif".format((len(parameters)-1)))
    
    #run elastix: 
    try:                
        if verbose: print ("Running Elastix, this can take some time....\n")
        sp.call(e_params)
        if verbose: print("Past Elastix Commandline Call")
    except RuntimeError as e:
        print("\n***RUNTIME ERROR***: {} Elastix has failed. Most likely the two images are too dissimiliar.\n".format(e.message))
        pass      
    if os.path.exists(ElastixResultFile) == True:    
        if verbose: print("Elastix Registration Successfully Completed\n")
    #check to see if it was MHD instead
    elif os.path.exists(os.path.join(out, "result.{}.mhd".format((len(parameters)-1)))) == True:    
        ElastixResultFile = os.path.join(out, "result.{}.mhd".format((len(parameters)-1)))
        if verbose: print("Elastix Registration Successfully Completed\n")
    else:
        print ("\n***ERROR***Cannot find elastix result file, try changing parameter files\n: {}".format(ElastixResultFile))
        return

    return ElastixResultFile, TransformParameterFile

if __name__ == "__main__":
    
    #setup
    print(sys.argv)
    stepid = int(sys.argv[1]) #0 or 1
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    #paths
    outdr = "/jukebox/scratch/zmd/registration_accuracy_iters/registration"
    src = "/jukebox/scratch/zmd/registration_accuracy_iters/volumes"
    vols = [os.path.join(src, xx) for xx in os.listdir(src)]
    
    #pick reg vol based on step ID (will only need 2) 
    start = time.time()
    vol = vols[stepid]
    fx = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
    out = os.path.join(outdr, "iter"+str(jobid).zfill(3))
    if not os.path.exists(out): os.mkdir(out)
    
    mv = vol
    
    params = ["/jukebox/wang/zahra/python/lightsheet_py3/parameterfolder/Order1_Par0000affine.txt", 
              "/jukebox/wang/zahra/python/lightsheet_py3/parameterfolder/Order2_Par0000bspline.txt"]
    
    elastix_command_line_call(fx, mv, out, params, fx_mask=False, verbose=False)
