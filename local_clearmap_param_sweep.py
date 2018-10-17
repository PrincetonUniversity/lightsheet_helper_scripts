#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 12:31:42 2018

@author: wanglab
"""

#Modifications of Tom's clearmap cluster code to parameter sweep an unprocessed brain locally

import os, sys, shutil
from xvfbwrapper import Xvfb; vdisplay = Xvfb(); vdisplay.start()
from ClearMap.cluster.preprocessing import updateparams, arrayjob, process_planes_completion_checker
from ClearMap.cluster.directorydeterminer import directorydeterminer
from ClearMap.cluster.par_tools import resampling_operations, celldetection_operations
from run_parameter_sweep import sweep_parameters_cluster
#%%

systemdirectory=directorydeterminer()
###set paths to data
###inputdictionary stucture: key=pathtodata value=list['xx', '##'] where xx=regch, injch, or cellch and ##=two digit channel number
#'regch' = channel to be used for registration, assumption is all other channels are signal
#'cellch' = channel(s) to apply cell detection
#'injch' = channels(s) to quantify injection site
#'##' = when taking a multi channel scan following regexpression, the channel corresponding to the reg/cell/inj channel. I.e. name_of_scan_channel00_Z#### then use '00'
#e.g.: inputdictionary={path_1: [['regch', '00']], path_2: [['cellch', '00'], ['injch', '01']]} ###create this dictionary variable BEFORE params
inputdictionary={
os.path.join(systemdirectory, 'LightSheetTransfer/Jess/cfos/180927_dadult_pc_crusi_7_053118_488_647_1d3x_017na_1hfds_z10um_150msec_15-35-08'): [['regch', '00'],['cellch', '01']]}

####Required inputs

######################################################################################################
#NOTE: edit clearmap/parameter_file.py for cell some detection parameters, everything else is handled below
######################################################################################################

params={
'inputdictionary': inputdictionary, #don't need to touch
'outputdirectory': os.path.join(systemdirectory, 'LightSheetTransfer/test/dadult_pc_crusi_7'),
'resample' : False, #False/None, float(e.g: 0.4), amount to resize by: >1 means increase size, <1 means decrease
'xyz_scale': (5.0, 5.0, 10), #micron/pixel; 1.3xobjective w/ 1xzoom 5um/pixel; 4x objective = 1.63um/pixel
'tiling_overlap': 0.00, #percent overlap taken during tiling
'AtlasFile' : os.path.join(systemdirectory, 'LightSheetTransfer/atlas/201810_jess_cfos/atlas_similarity_to_auto_masked.tif'), ###it is assumed that input image will be a horizontal scan with anterior being 'up'; USE .TIF!!!!
'annotationfile' :   os.path.join(systemdirectory, 'LightSheetTransfer/atlas/201810_jess_cfos/ann_similarity_to_auto_masked.tif'), ###path to annotation file for structures
'blendtype' : 'sigmoidal', #False/None, 'linear', or 'sigmoidal' blending between tiles, usually sigmoidal; False or None for images where blending would be detrimental;
'intensitycorrection' : False, #True = calculate mean intensity of overlap between tiles shift higher of two towards lower - useful for images where relative intensity is not important (i.e. tracing=True, cFOS=False)
'rawdata' : True, # set to true if raw data is taken from scope and images need to be flattened; functionality for rawdata =False has not been tested**
'FinalOrientation': (3, 2, 1), #Orientation: 1,2,3 means the same orientation as the reference and atlas files; #Flip axis with - sign (eg. (-1,2,3) flips x). 3D Rotate by swapping numbers. (eg. (2,1,3) swaps x and y); USE (3,2,1) for DVhorizotnal to sagittal. NOTE (TP): -3 seems to mess up the function and cannot seem to figure out why. do not use.
'slurmjobfactor': 650, #number of array iterations per arrayjob since max job array on SPOCK is 1000
'removeBackgroundParameter_size': (5,5), #Remove the background with morphological opening (optimised for spherical objects), e.g. (7,7)
'findExtendedMaximaParameter_hmax': None, # (float or None)     h parameter (for instance 20) for the initial h-Max transform, if None, do not perform a h-max transform
'findExtendedMaximaParameter_size': 5, # size in pixels (x,y) for the structure element of the morphological opening
'findExtendedMaximaParameter_threshold': 0, # (float or None)     include only maxima larger than a threshold, if None keep all local maxima
'findIntensityParameter_method': 'Max', # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
'findIntensityParameter_size': (3,3,3), # (tuple)             size of the search box on which to perform the *method*
'detectCellShapeParameter_threshold': 105 # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated
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

#to run locally:
#from ClearMap.cluster.process_local import run_brain_locally
#run_brain_locally(steps = [4,5,6], **params)


#run scipt portions        arrayjob(jobid, cores=5, compression=1, **params) #process zslice numbers equal to slurmjobfactor*jobid thru (jobid+1)*slurmjobfactor

if __name__ == '__main__':

    #get job id from SLURM
    print sys.argv
    print os.environ["SLURM_ARRAY_TASK_ID"]
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"]) #int(sys.argv[2])
    stepid = int(sys.argv[1])

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
    arrayjob(0, cores=5, compression=1, **params) #process zslice numbers equal to slurmjobfactor*jobid thru (jobid+1)*slurmjobfactor
    arrayjob(1, cores=5, compression=1, **params)
        
    #######################STEP 2 #######################
    #####################################################
    ###check to make sure all step 1 jobs completed properly
    process_planes_completion_checker(**params)
    #clearmap: load the parameters:
    resampling_operations(0, **params)
    resampling_operations(1, **params)


#%%
##add way to easily test cell detection parameters. Just select a 'jobid': e.g. jobid = 22
    #parameter sweep cell detection parameters. NOTE read all of functions description before using. VERY CPU intensive
    for jobid in range(5):
        try:
            #parameter sweep cell detection parameters. NOTE read all of functions description before using. VERY CPU intensive
            sweep_parameters_cluster(jobid, cleanup = False, **params)
        except:
            try:
                ###make parameter dictionary and pickle file:
                updateparams(os.getcwd(), **params) # e.g. single job assuming directory_determiner function has been properly set
                #copy folder into output for records
                if not os.path.exists(os.path.join(params['outputdirectory'], 'clearmap_cluster')): shutil.copytree(os.getcwd(), os.path.join(params['outputdirectory'], 'clearmap_cluster'), ignore=shutil.ignore_patterns('^.git')) #copy run folder into output to save run info
                sweep_parameters_cluster(jobid, cleanup = False , **params)
            except Exception, e:
                print('Jobid {}, Error given {}'.format(jobid, e))
