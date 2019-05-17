#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 18:05:44 2019

@author: wanglab
"""

import os, sys, shutil
from xvfbwrapper import Xvfb; vdisplay = Xvfb(); vdisplay.start()
sys.path.append("/jukebox/wang/zahra/lightsheet_copy")
from tools.imageprocessing import preprocessing
from tools.utils.directorydeterminer import directorydeterminer

systemdirectory=directorydeterminer()

###set paths to data
###inputdictionary stucture: key=pathtodata value=list['xx', '##'] where xx=regch, injch, or cellch and ##=two digit channel number
#'regch' = channel to be used for registration, assumption is all other channels are signal
#'cellch' = channel(s) to apply cell detection
#'injch' = channels(s) to quantify injection site
#e.g.: inputdictionary={path_1: [['regch', '00']], path_2: [['cellch', '00'], ['injch', '01']]} ###create this dictionary variable BEFORE params
inputdictionary={
os.path.join(systemdirectory, 'LightSheetTransfer/kelly/20190415_cfos_for_microscope_tests/190505_m57207_dems_cfos_20190320_4x_647_008na_1hfds_z2um_200msec_10povlp_12-38-25'): [['cellch', '00']]}
####Required inputs

params={
'labeltype': 'h129', #'h129', 'prv', 'cfos'
'objectdetection': 'convnet', # 'edgedetection', 'convnet', 'clearmap', 'all'; clearmap setting uses their SpotDetection method
'systemdirectory':  systemdirectory, #don't need to touch
'inputdictionary': inputdictionary, #don't need to touch
'outputdirectory': os.path.join(systemdirectory, 'LightSheetTransfer/kelly/20190415_cfos_for_microscope_tests/m57207_dems_cfos_20190320_4x_stitched'),
'xyz_scale': (1.63, 1.63, 2), #(5.0,5.0,3), #micron/pixel: 5.0um/pix for 1.3x; 1.63um/pix for 4x
'tiling_overlap': 0.10, #percent overlap taken during tiling
'stitchingmethod': 'terastitcher', #'terastitcher', blending see below for details
'AtlasFile' : os.path.join(systemdirectory, 'wang/pisano/Python/atlas/sagittal_atlas_20um_iso.tif'),
'annotationfile' : os.path.join(systemdirectory, 'wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif'), ###path to annotation file for structures
'blendtype' : 'sigmoidal', #False/None, 'linear', or 'sigmoidal' blending between tiles, usually sigmoidal; False or None for images where blending would be detrimental
'intensitycorrection' : True, #True = calculate mean intensity of overlap between tiles shift higher of two towards lower - useful for images where relative intensity is not important (i.e. tracing=True, cFOS=False)
'resizefactor': 6, ##in x and y #normally set to 5 for 4x objective, 3 for 1.3x obj
'rawdata' : True, # set to true if raw data is taken from scope and images need to be flattened; functionality for rawdata =False has not been tested**
'finalorientation' :  ('2','1','0'), #Used to account for different orientation between brain and atlas. Assumes XYZ ('0','1','2) orientation. Pass strings NOT ints. '-0' = reverse the order of the xaxis. For better description see docstring from tools.imageprocessing.orientation import fix_orientation; ('2','1','0') for horizontal to sagittal, Order of operations is reversing of axes BEFORE swapping axes.
'slurmjobfactor': 20 #number of array iterations per arrayjob since max job array on SPOCK is 1000
}
###make parameter dictionary and pickle file:
#preprocessing.generateparamdict(os.getcwd(), **params) # e.g. single job assuming directory_determiner function has been properly set
##preprocessing.updateparams('/home/wanglab/wang/pisano/Python/lightsheet', svnm = 'param_dict_local.p', **params) # make a local copy
#if not os.path.exists(os.path.join(params['outputdirectory'], 'lightsheet')): shutil.copytree(os.getcwd(), os.path.join(params['outputdirectory'], 'lightsheet'), ignore=shutil.ignore_patterns(*('.pyc','CVS','.git','tmp','.svn', 'TeraStitcher-Qt4-standalone-1.10.11-Linux'))) #copy run folder into output to save run info
##os.system("rsync -av --exclude='.git/' ....)#
#######################STEP 1 #######################
#Stitch and preprocess each z plane
#####################################################
#Stitch using Terastitcher "smart stitching"
from tools.imageprocessing.stitch import terastitcher_from_params
terastitcher_from_params(**params)

preprocessing.tiffcombiner(0, **params)
