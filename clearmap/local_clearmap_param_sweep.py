#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 12:31:42 2018

@author: wanglab
"""

#Modifications of Tom's clearmap cluster code to parameter sweep an unprocessed brain locally

import os, sys, shutil, tifffile, numpy as np
from skimage.exposure import rescale_intensity
from itertools import product
from ClearMap.cluster.preprocessing import updateparams, arrayjob, process_planes_completion_checker
from ClearMap.cluster.directorydeterminer import directorydeterminer
from ClearMap.cluster.par_tools import resampling_operations, celldetection_operations
from ClearMap.cluster.preprocessing import makedir, listdirfull, removedir
from ClearMap.cluster.utils import load_kwargs


systemdirectory=directorydeterminer()
###set paths to data
###inputdictionary stucture: key=pathtodata value=list['xx', '##'] where xx=regch, injch, or cellch and ##=two digit channel number
#'regch' = channel to be used for registration, assumption is all other channels are signal
#'cellch' = channel(s) to apply cell detection
#'injch' = channels(s) to quantify injection site
#'##' = when taking a multi channel scan following regexpression, the channel corresponding to the reg/cell/inj channel. I.e. name_of_scan_channel00_Z#### then use '00'
#e.g.: inputdictionary={path_1: [['regch', '00']], path_2: [['cellch', '00'], ['injch', '01']]} ###create this dictionary variable BEFORE params
inputdictionary={
os.path.join(systemdirectory, 'LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/raw_data/190123_buffer_retroorb_488_647_017na_1hfsrs_z10um_75ms_14-57-55'): [['regch', '00'],['cellch', '01']]
 }

####Required inputs

######################################################################################################
#NOTE: edit clearmap/parameter_file.py for cell some detection parameters, everything else is handled below
######################################################################################################

params={
'inputdictionary': inputdictionary, #don't need to touch
'outputdirectory': os.path.join(systemdirectory, 'LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/parameter_sweep/buffer'),
'resample' : False, #False/None, float(e.g: 0.4), amount to resize by: >1 means increase size, <1 means decrease
'xyz_scale': (5.0, 5.0, 10), #micron/pixel; 1.3xobjective w/ 1xzoom 5um/pixel; 4x objective = 1.63um/pixel
'tiling_overlap': 0.00, #percent overlap taken during tiling
'AtlasFile' : os.path.join(systemdirectory, 'LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/average_template_25_sagittal_forDVscans_z_thru_240.tif'), ###it is assumed that input image will be a horizontal scan with anterior being 'up'; USE .TIF!!!!
'annotationfile' :   os.path.join(systemdirectory, 'LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/atlas/annotation_25_ccf2015_forDVscans_z_thru_240.nrrd'), ###path to annotation file for structures
'blendtype' : 'sigmoidal', #False/None, 'linear', or 'sigmoidal' blending between tiles, usually sigmoidal; False or None for images where blending would be detrimental;
'intensitycorrection' : False, #True = calculate mean intensity of overlap between tiles shift higher of two towards lower - useful for images where relative intensity is not important (i.e. tracing=True, cFOS=False)
'rawdata' : True, # set to true if raw data is taken from scope and images need to be flattened; functionality for rawdata =False has not been tested**
'FinalOrientation': (3, 2, 1), #Orientation: 1,2,3 means the same orientation as the reference and atlas files; #Flip axis with - sign (eg. (-1,2,3) flips x). 3D Rotate by swapping numbers. (eg. (2,1,3) swaps x and y); USE (3,2,1) for DVhorizotnal to sagittal. NOTE (TP): -3 seems to mess up the function and cannot seem to figure out why. do not use.
'slurmjobfactor': 50, #number of array iterations per arrayjob since max job array on SPOCK is 1000
'removeBackgroundParameter_size': (5,5), #Remove the background with morphological opening (optimised for spherical objects), e.g. (7,7)
'findExtendedMaximaParameter_hmax': None, # (float or None)     h parameter (for instance 20) for the initial h-Max transform, if None, do not perform a h-max transform
'findExtendedMaximaParameter_size': 5, # size in pixels (x,y) for the structure element of the morphological opening
'findExtendedMaximaParameter_threshold': 0, # (float or None)     include only maxima larger than a threshold, if None keep all local maxima
'findIntensityParameter_method': 'Max', # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
'findIntensityParameter_size': (3,3,3), # (tuple)             size of the search box on which to perform the *method*
'detectCellShapeParameter_threshold': 125 # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated
}

#####################################################################################################################################################
##################################################optional arguments for params######################################################################
#####################################################################################################################################################
#'regexpression':  r'(.*)(?P<y>\d{2})(.*)(?P<x>\d{2})(.*C+)(?P<ch>[0-9]{1,2})(.*Z+)(?P<z>[0-9]{1,4})(.ome.tif)', ###lavision preprocessed data
#'regexpression':  r'(.*)(.*C+)(?P<ch>[0-9]{1,2})(.*Z+)(?P<z>[0-9]{1,4})(.ome.tif)', lavision NONTILING**
#regexpression: 'r'(.*)(.*C+)(?P<ch>[0-9]{1,2})(.*Z+)(?P<z>[0-9]{1,4})(.ome.tif)', lavision Channels and Z
#'parameterfolder' : os.path.join(systemdirectory, 'wang/pisano/Python/lightsheet/parameterfolder'), ##  * folder consisting of elastix parameter files with prefixes "Order<#>_" to specify application order
#'removeBackgroundParameter_size': (7,7), #Remove the background with morphological opening (optimised for spherical objects), e.g. (7,7)
#'findExtendedMaximaParameter_hmax': None, # (float or None)     h parameter (for instance 20) for the initial h-Max transform, if None, do not perform a h-max transform
#'findExtendedMaximaParameter_size': 5 # size in pixels (x,y) for the structure element of the morphological opening
#'findExtendedMaximaParameter_threshold': 0, # (float or None)     include only maxima larger than a threshold, if None keep all local maxima
#'findIntensityParameter_method': 'Max', # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
#'findIntensityParameter_size': (3,3,3), # (tuple)             size of the search box on which to perform the *method*
#'detectCellShapeParameter_threshold': 500 # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated
#####################################################################################################################################################
#####################################################################################################################################################
#####################################################################################################################################################

#sweep parameters copy & modifications - run before running above cell                
def sweep_parameters_cluster(jobid, optimization_chunk=7, pth=False, rescale=False, cleanup=True, **kwargs):
    '''Function to sweep parameters
    
    final outputs will be saved in outputdirectory/parameter_sweep
    second copy will be saved in outputdirectory/parameter_sweep_jobid if cleanup=False

    Inputs:
        ----------------
        jobid: chunk of tissue to run (usually int between 20-30)
        #pth (optional): if pth to output folder after running package, function will load the param file automatically
        rescale (optional): str of dtype to rescale to. E.g.: 'uint8'
        cleanup = T/F removes subfolders after
        optimization_chunk = this was the old "jobid" in this case it is the chunk of volume to look at
        kwargs (if not pth): 'params' from run_clearmap_cluster.py
    '''

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
#    dCSP_threshold_r = range(325,450,25) #<-- IMPORTANT TO SWEEP
#    
    #second cleanup=False
    rBP_size_r = [3] #zmd commented out
    dCSP_threshold_r = [440, 450, 460, 470, 480, 490, 500]
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    
    
    # calculate number of iterations
    tick = 0
    for rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold in product(rBP_size_r, fEMP_hmax_r, fEMP_size_r, fEMP_threshold_r, fIP_method_r, fIP_size_r, dCSP_threshold_r):
        tick +=1

    sys.stdout.write('\n\nNumber of iterations is {}:'.format(tick))
    
    #if pth is set - zmd added
    if pth:
        kwargs = load_kwargs(pth)
        
    #make folder for final output:
    opt = kwargs['outputdirectory']; makedir(opt)
    out = opt+'/parameter_sweep'; makedir(out)
    out0 = opt+'/parameter_sweep_jobid_{}'.format(str(jobid).zfill(4)); makedir(out0)

    ntick = 0
    rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold=[(rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold) for rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold in product(rBP_size_r, fEMP_hmax_r, fEMP_size_r, fEMP_threshold_r, fIP_method_r, fIP_size_r, dCSP_threshold_r)][jobid]

    #zmd modified
    pth = out0+'/parametersweep_rBP_size{}_fEMP_hmax{}_fEMP_size{}_fEMP_threshold{}_fIP_method{}_fIP_size{}_dCSP_threshold{}.tif'.format(rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold)

    if not os.path.exists(pth):

        try:

            #set params for sweep
            kwargs['removeBackgroundParameter_size'] = (rBP_size, rBP_size) #Remove the background with morphological opening (optimised for spherical objects), e.g. (7,7)
            kwargs['findExtendedMaximaParameter_hmax'] = fEMP_hmax # (float or None)     h parameter (for instance 20) for the initial h-Max transform, if None, do not perform a h-max transform
            kwargs['findExtendedMaximaParameter_size'] = fEMP_size # size in pixels (x,y) for the structure element of the morphological opening
            kwargs['findExtendedMaximaParameter_threshold'] = fEMP_threshold # (float or None)     include only maxima larger than a threshold, if None keep all local maxima
            kwargs['findIntensityParameter_method'] =  fIP_method # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
            kwargs['findIntensityParameter_size'] = (fIP_size,fIP_size,fIP_size) # (tuple)             size of the search box on which to perform the *method*
            kwargs['detectCellShapeParameter_threshold'] = dCSP_threshold # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated

            #tmp
            import cPickle as pickle
            from ClearMap.cluster.utils import load_kwargs
            nkwargs = load_kwargs(kwargs['outputdirectory'])
            kwargs['outputdirectory'] = out0
            nkwargs.update(kwargs)
            pckloc=out0+'/param_dict.p'; pckfl=open(pckloc, 'wb'); pickle.dump(nkwargs, pckfl); pckfl.close()

            #run cell detection
            ntick+=1
            sys.stdout.write('\n\n\n           *****Iteration {} of {}*****\n\n\n'.format(ntick, tick))
            sys.stdout.write('    Iteration parameters: {}     {}     {}     {}     {}     {}     {}'.format(kwargs['removeBackgroundParameter_size'], kwargs['findExtendedMaximaParameter_hmax'], kwargs['findExtendedMaximaParameter_size'], kwargs['findExtendedMaximaParameter_threshold'],         kwargs['findIntensityParameter_method'],         kwargs['findIntensityParameter_size'],        kwargs['detectCellShapeParameter_threshold']))
            celldetection_operations(optimization_chunk, testing = True, **kwargs)

            #list, load, and maxip
            if ntick == 1: raw = [xx for xx in listdirfull(out0+'/optimization/raw') if '~' not in xx and '.db' not in xx]; raw.sort(); raw_im = np.squeeze(tifffile.imread(raw)); raw_mx = np.max(raw_im, axis = 0)
            bkg = [xx for xx in listdirfull(out0+'/optimization/background') if '~' not in xx and 'Thumbs.db' not in xx]; bkg.sort(); bkg_im = tifffile.imread(bkg); bkg_mx = np.max(bkg_im, axis = 0)
            cell = [xx for xx in listdirfull(out0+'/optimization/cell') if '~' not in xx and '.db' not in xx]; cell.sort(); cell_im = tifffile.imread(cell); cell_mx = np.max(cell_im, axis = 0)

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
            npth = out+'/jobid_{}_parametersweep_rBP_size{}_fEMP_hmax{}_fEMP_size{}_fEMP_threshold{}_fIP_method{}_fIP_size{}_dCSP_threshold{}.tif'.format(str(jobid).zfill(4), rBP_size, fEMP_hmax, fEMP_size, fEMP_threshold, fIP_method, fIP_size, dCSP_threshold)
            tifffile.imsave(npth, bigim, compress = 1)
            

        except Exception, e:
            print ('Error on: {}\n\nerror={}'.format(pth,e))
            im = np.zeros((10,10,10))
            tifffile.imsave(pth, im, compress = 1)
            with open(os.path.join(out, 'errored_files.txt'), 'a') as fl:
                fl.write('\n\n{}\n{}\n'.format(pth, kwargs))
                fl.close

    return

#%%
#run scipt portions        
#this particular cell only stitches (by blending) and resamples the data

if __name__ == '__main__':
    
    #######################STEP 0 #######################
    #####################################################
    ###make parameter dictionary and pickle file:
    updateparams(os.getcwd(), **params) # e.g. single job assuming directory_determiner function has been properly set
    #copy folder into output for records
    if not os.path.exists(os.path.join(params['outputdirectory'], 'clearmap_cluster')): shutil.copytree(os.getcwd(), os.path.join(params['outputdirectory'], 'clearmap_cluster'), ignore=shutil.ignore_patterns('^.git')) #copy run folder into output to save run info

    #re-run update params like you would on the cluster
    updateparams(os.path.join(params['outputdirectory'], 'clearmap_cluster'), **params)

    #######################STEP 1 #######################
    #####################################################
    ###stitch, resample, and save files
    arrayjob(0, cores=8, compression=1, **params) #process zslice numbers equal to slurmjobfactor*jobid thru (jobid+1)*slurmjobfactor
    arrayjob(1, cores=8, compression=1, **params)
        
    #######################STEP 2 #######################
    #####################################################
    ###check to make sure all step 1 jobs completed properly
    process_planes_completion_checker(**params)
    #clearmap: load the parameters:
    resampling_operations(0, **params)
    resampling_operations(1, **params)


    ##add way to easily test cell detection parameters. Just select a 'jobid': e.g. jobid = 22
    #parameter sweep cell detection parameters. NOTE read all of functions description before using. VERY CPU intensive
    #for first pass at cell detection
    for jobid in range(7): #to find the range value, run lines 152 - 173 in this script, part of sweep_parameters_cluster func
        
        sweep_parameters_cluster(jobid, **params)
            