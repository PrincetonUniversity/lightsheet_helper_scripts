#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 12:43:11 2019

@author: wanglab
"""
import os, sys, numpy as np, subprocess as sp
from skimage.external import tifffile
import numpy as np
from scipy.ndimage import zoom

def elastix_command_line_call(fx, mv, out, parameters, fx_mask=False, verbose=False):
    '''Wrapper Function to call elastix using the commandline, this can be time consuming
    
    Inputs
    -------------------
    fx = fixed path (usually Atlas for 'normal' noninverse transforms)
    mv = moving path (usually volume to register for 'normal' noninverse transforms)
    out = folder to save file
    parameters = list of paths to parameter files IN ORDER THEY SHOULD BE APPLIED
    fx_mask= (optional) mask path if desired
    
    Outputs
    --------------
    ElastixResultFile = '.tif' or '.mhd' result file
    TransformParameterFile = file storing transform parameters
    
    '''
    e_params=['elastix', '-f', fx, '-m', mv, '-out', out]
    if fx_mask: e_params=['elastix', '-f', fx, '-m', mv, '-fMask', fx_mask, '-out', out]
    
    ###adding elastix parameter files to command line call
    for x in range(len(parameters)):
        e_params.append('-p')
        e_params.append(parameters[x])
    
    #set paths
    TransformParameterFile = os.path.join(out, 'TransformParameters.{}.txt'.format((len(parameters)-1)))
    ElastixResultFile = os.path.join(out, 'result.{}.tif'.format((len(parameters)-1)))
    
    #run elastix: 
    try:                
        if verbose: print ('Running Elastix, this can take some time....\n')
        sp.call(e_params)
        if verbose: print('Past Elastix Commandline Call')
    except RuntimeError as e:
        print('\n***RUNTIME ERROR***: {} Elastix has failed. Most likely the two images are too dissimiliar.\n'.format(e.message))
        pass      
    if os.path.exists(ElastixResultFile) == True:    
        if verbose: print('Elastix Registration Successfully Completed\n')
    #check to see if it was MHD instead
    elif os.path.exists(os.path.join(out, 'result.{}.mhd'.format((len(parameters)-1)))) == True:    
        ElastixResultFile = os.path.join(out, 'result.{}.mhd'.format((len(parameters)-1)))
        if verbose: print('Elastix Registration Successfully Completed\n')
    else:
        print ('\n***ERROR***Cannot find elastix result file, try changing parameter files\n: {}'.format(ElastixResultFile))
        return

        
    return ElastixResultFile, TransformParameterFile

#%%
#using registration channel of ai148 with a lot of PC staining
#only downsized volume, not resampled
dwnsz = "/jukebox/wang/mkislin/lightsheet_brains/201903_cntnap2_tsc1_ai148/ai148_47018_ii/20190130_mk_ai148_47018_ii_1d3x_488_647_008na_1hfds_z10um_200msec_resized_ch00.tif"

##resampling to pma atlas coordinates
#dwnsz = tifffile.imread(dwnsz)
#factor = (702/dwnsz.shape[0], 832/dwnsz.shape[1], 457/dwnsz.shape[2]) #got these coordinates from a tracing brain registered to pma
#mv = zoom(dwnsz, factor, order = 3)
##save out
#tifffile.imsave("/jukebox/wang/zahra/modeling/cerebellar_cortex_segmentation/ai148_47018_ii_resampled.tif", mv)
mv = "/jukebox/wang/zahra/modeling/cerebellar_cortex_segmentation/ai148_47018_ii_resampled.tif"
#pma atlas
fx = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"

#setup destination
dst = "/jukebox/wang/zahra/modeling/cerebellar_cortex_segmentation/ai148"

#registration params for mouse
params = ["/jukebox/wang/zahra/lightsheet_copy/parameterfolder/Order1_Par0000affine.txt", 
          "/jukebox/wang/zahra/lightsheet_copy/parameterfolder/Order2_Par0000bspline.txt"]

#register
elastix_command_line_call(fx = fx, mv = mv, out = dst, parameters = params, fx_mask=False)

#now transform cell channel to registration channel?

