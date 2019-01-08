#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 11:27:50 2018

@author: wanglab
"""

import os, shutil, re
from xvfbwrapper import Xvfb; vdisplay = Xvfb(); vdisplay.start()
from tools.imageprocessing import preprocessing
from tools.registration.register import elastix_wrapper

#mimics run_tracing file with cerebellum-specific parameters to register locally.
#does not perform elastix inverse transform; registration primarily for using analyze_injection function on injection site data.

#cerebellums to process
inputs = [
        '/jukebox/LightSheetTransfer/Jess/181221_pcdev_lob6_2_cb_488_75p_647_1d3x_026na_1hfds_z10um_250msec_15-14-26'
        ]

#%%
for brain in inputs: #for loop that processes each brain one by one
    
    inputdictionary = { brain: [['regch', '00'],['injch', '01']] } #makes input directory 
    
    brainname = re.search('(?<=_)(\w+)(?=_488)', brain) #sets basename of output directory based on clearing label
    
    params={
    'labeltype': 'CTB555', #'h129', 'prv', 'cfos'
    'objectdetection': 'edgedetection', # 'edgedetection', 'convnet', 'clearmap', 'all'; clearmap setting uses their SpotDetection method
    'systemdirectory':  '/jukebox/', #changed for local processing by zmd - confuses file paths less
    'inputdictionary': inputdictionary, #don't need to touch
    'outputdirectory': '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/'+brainname.group(0),
    'xyz_scale': (5,5,10), #(5.0,5.0,3), #micron/pixel: 5.0um/pix for 1.3x; 1.63um/pix for 4x; The third number, Z, is the size of the z-step
    'tiling_overlap': 0.00, #percent overlap taken during tiling
    'stitchingmethod': 'blending', #'terastitcher', blending see below for details
    'AtlasFile' : '/jukebox/LightSheetTransfer/atlas/cb_sagittal_atlas_20um_iso.tif',
    'annotationfile' : '/jukebox/LightSheetTransfer/atlas/cb_annotation_sagittal_atlas_20um_iso.tif', ###path to annotation file for structures
    'blendtype' : 'sigmoidal', #False/None, 'linear', or 'sigmoidal' blending between tiles, usually sigmoidal; False or None for images where blending would be detrimental
    'intensitycorrection' : True, #True = calculate mean intensity of overlap between tiles shift higher of two towards lower - useful for images where relative intensity is not important (i.e. tracing=True, cFOS=False)
    'resizefactor': 3, ##in x and y #normally set to 5 for 4x objective, 3 for 1.3x obj
    'rawdata' : True, # set to true if raw data is taken from scope and images need to be flattened; functionality for rawdata =False has not been tested**
    'finalorientation' :  ('2','1','0'), #Used to account for different orientation between brain and atlas. Assumes XYZ ('0','1','2) orientation. Pass strings NOT ints. '-0' = reverse the order of the xaxis. For better description see docstring from tools.imageprocessing.orientation import fix_orientation; ('2','1','0') for horizontal to sagittal, Order of operations is reversing of axes BEFORE swapping axes.
    'slurmjobfactor': 700, #number of array iterations per arrayjob since max job array on SPOCK is 1000
    'transfertype': 'copy', #to protect original data
    'secondary_registration': False #to prevent morphing of signal channel image with unnecessary registration
    } 

    for stepid in [0,1,11,2,3]: #iterates through only relevant steps
        
        #######################STEP 0 #######################
        #Make parameter dictionary and setup destination
        #####################################################
        if stepid == 0:
            ###make parameter dictionary and pickle file:
            preprocessing.generateparamdict(os.getcwd(), **params) # e.g. single job assuming directory_determiner function has been properly set

            if not os.path.exists(os.path.join(params['outputdirectory'], 'lightsheet')): shutil.copytree(os.getcwd(), os.path.join(params['outputdirectory'], 'lightsheet'), ignore=shutil.ignore_patterns(*('.pyc','CVS','.git','tmp','.svn', 'TeraStitcher-Qt4-standalone-1.10.11-Linux'))) #copy run folder into output to save run info
            
            #then, re-run step 0 like you would do on the cluster to set up the proper param_dict file based on outputdirectory
            preprocessing.generateparamdict(os.path.join(params['outputdirectory'], 'lightsheet'), **params)  #re-run this as if you would do on the cluster
                    
        #######################STEP 1 #######################
        #Stitch and preprocess each z plane
        #####################################################
        elif stepid == 1:
            ###stitch based on percent overlap only ('dumb stitching'), and save files; showcelldetection=True: save out cells contours ovelaid on images
            preprocessing.arrayjob(0, cores = 1, compression = 1, **params) #process zslice numbers equal to slurmjobfactor*jobid thru (jobid+1)*slurmjobfactor
                
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
           elastix_wrapper(0, cores = 12, **params) #run elastix    