#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 11:27:50 2018

@author: wanglab
"""

import os, sys, shutil
from xvfbwrapper import Xvfb; vdisplay = Xvfb(); vdisplay.start()
from tools.imageprocessing import preprocessing
from tools.registration.register import elastix_wrapper
from tools.utils.directorydeterminer import directorydeterminer
from tools.utils.io import load_kwargs

#cerebellums to process
inputs = [
        '/jukebox/LightSheetTransfer/Jess/Lawrence/180706_lawrence_an10_crus_iDisco_488_647_025na_1hfds_z10um_250msec_18-27-43',
        '/jukebox/LightSheetTransfer/Jess/Lawrence/180706_lawrence_an11_crus_iDisco_488_647_025na_1hfds_z10um_250msec_18-50-49',
        '/jukebox/LightSheetTransfer/Jess/Lawrence/180706_lawrence_an12_crus_iDisco_488_647_025na_1hfds_z10um_250msec_19-13-20'
        ]
        #'/jukebox/LightSheetTransfer/Jess/Lawrence/180706_lawrence_an7_crus_iDisco_488_647_025na_1hfds_z10um_250msec_17-05-32',
        #'/jukebox/LightSheetTransfer/Jess/christina/lobVI_DREADD/180706_christina_an21_iDisco_488_647_025na_1hfds_z10um_250msec_14-59-47',
        #'/jukebox/LightSheetTransfer/Jess/christina/crusI_DREADD/180828_christina_an2_iDisco_488_647_026na_1hfds_z10um_250msec_13-34-54'
        


for brain in inputs: #for loop that processes each brain one by one
    inputdictionary = { brain: [['regch', '00'],['injch', '01']] } #makes input directory 
    
    brainname = brain[49:108]
    #if brain[63:118].endswith('_'): #FIXME: use regex to do this!!
    #    brainname = brain[63:117]
    #else: 
    #    brainname = brain[63:118]
    
    params={
    'labeltype': 'CTB555', #'h129', 'prv', 'cfos'
    'objectdetection': 'edgedetection', # 'edgedetection', 'convnet', 'clearmap', 'all'; clearmap setting uses their SpotDetection method
    'systemdirectory':  '/jukebox/', #don't need to touch
    'inputdictionary': inputdictionary, #don't need to touch
    'outputdirectory': '/jukebox/wang/Jess/lightsheet_output/lawrence/'+brainname,
    'xyz_scale': (5,5,10), #(5.0,5.0,3), #micron/pixel: 5.0um/pix for 1.3x; 1.63um/pix for 4x; The third number, Z, is the size of the z-step
    'tiling_overlap': 0.00, #percent overlap taken during tiling
    'stitchingmethod': 'blending', #'terastitcher', blending see below for details
    'AtlasFile' : '/jukebox/LightSheetTransfer/atlas/cb_sagittal_atlas_20um_iso.tif',
    'annotationfile' : '/jukebox/LightSheetTransfer/atlas/cb_annotation_sagittal_atlas_20um_iso.tif', ###path to annotation file for structures
    'blendtype' : 'sigmoidal', #False/None, 'linear', or 'sigmoidal' blending between tiles, usually sigmoidal; False or None for images where blending would be detrimental
    'intensitycorrection' : True, #True = calculate mean intensity of overlap between tiles shift higher of two towards lower - useful for images where relative intensity is not important (i.e. tracing=True, cFOS=False)
    'resizefactor': 3, ##in x and y #normally set to 5 for 4x objective, 3 for 1.3x obj
    'rawdata' : True, # set to true if raw data is taken from scope and images need to be flattened; functionality for rawdata =False has not been tested**
    'finalorientation' :  ('2','-1','0'), #Used to account for different orientation between brain and atlas. Assumes XYZ ('0','1','2) orientation. Pass strings NOT ints. '-0' = reverse the order of the xaxis. For better description see docstring from tools.imageprocessing.orientation import fix_orientation; ('2','1','0') for horizontal to sagittal, Order of operations is reversing of axes BEFORE swapping axes.
    'slurmjobfactor': 600, #number of array iterations per arrayjob since max job array on SPOCK is 1000
    'transfertype': 'copy', #to protect original data
    'secondary_registration': False #to prevent morphing of signal channel image with unnecessary registration
    } 
    

    for stepid in [0,1,11,2,3]:
        
        #######################STEP 0 #######################
        #Make parameter dictionary and setup destination
        #####################################################
        if stepid == 0:
            ###make parameter dictionary and pickle file:
            preprocessing.generateparamdict(os.getcwd(), **params) # e.g. single job assuming directory_determiner function has been properly set

            if not os.path.exists(os.path.join(params['outputdirectory'], 'lightsheet')): shutil.copytree(os.getcwd(), os.path.join(params['outputdirectory'], 'lightsheet'), ignore=shutil.ignore_patterns(*('.pyc','CVS','.git','tmp','.svn', 'TeraStitcher-Qt4-standalone-1.10.11-Linux'))) #copy run folder into output to save run info
            
            #kwargs = load_kwargs(params['outputdirectory'])#re-load incorrect kwargs
                      
            preprocessing.generateparamdict(os.path.join(params['outputdirectory'], 'lightsheet'), **params)  #re-run this as if you would do on the cluster
                    
        #######################STEP 1 #######################
        #Stitch and preprocess each z plane
        #####################################################
        elif stepid == 1:
            ###stitch based on percent overlap only ('dumb stitching'), and save files; showcelldetection=True: save out cells contours ovelaid on images
            preprocessing.arrayjob(0, cores=1, compression=1, **params) #process zslice numbers equal to slurmjobfactor*jobid thru (jobid+1)*slurmjobfactor
                
        #######################STEP 1 check##################
        #Check to make sure jobs weren't lost - important if running on cluster
        #####################################################
        elif stepid == 11:
            ###check to make sure all step 1 jobs completed properly
            preprocessing.process_planes_completion_checker(**params)
        
        #######################STEP 2 #######################
        #Consolidate for Registration
        #####################################################
        elif stepid == 2:
            ###combine downsized ch and ch+cell files
            preprocessing.tiffcombiner(0, cores = 8, **params)
            preprocessing.tiffcombiner(1, cores = 8, **params)
            
        #######################STEP 3 #######################
        #####################################################
        elif stepid == 3:
           elastix_wrapper(0, cores = 6, **params) #run elastix    