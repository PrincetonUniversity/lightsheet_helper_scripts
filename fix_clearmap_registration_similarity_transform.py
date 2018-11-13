#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:14:44 2018

@author: wanglab
"""

import os
from ClearMap.cluster.utils import load_kwargs
from ClearMap.cluster.preprocessing import makedir
os.chdir('/jukebox/wang/zahra/lightsheet_copy')
from tools.registration.registration_using_similarity_mask import mask_similarity_transformed_atlas
from tools.utils.io import listdirfull, removedir, chunkit, writer, convert_to_mhd
import subprocess as sp

ann = '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif'
atlas = '/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif'
src = '/jukebox/wang/Jess/lightsheet_output/201810_cfos/processed/dadult_pc_crus1_1'

def elastix_command_line_call(fx, mv, out, parameters, fx_mask=False):
    '''Wrapper Function to call elastix using the commandline, this can be time consuming
    
    by @tpisano
    '''
    e_params=['elastix', '-f', fx, '-m', mv, '-out', out]
    if fx_mask: e_params=['elastix', '-f', fx, '-m', mv, '-fMask', fx_mask, '-out', out]
    
    ###adding elastix parameter files to command line call
    for x in range(len(parameters)):
        e_params.append('-p')
        e_params.append(parameters[x])
    writer(out,'Elastix Command:\n{}\n...'.format(e_params))    
    
    #set paths
    TransformParameterFile = os.path.join(out, 'TransformParameters.{}.txt'.format((len(parameters)-1)))
    ElastixResultFile = os.path.join(out, 'result.{}.tif'.format((len(parameters)-1)))
    
    #run elastix: 
    try:                
        print ('Running Elastix, this can take some time....\n')
        sp.call(e_params)#sp_call(e_params)#
        writer(out,'Past Elastix Commandline Call')
    except RuntimeError, e:
        writer(out,'\n***RUNTIME ERROR***: {} Elastix has failed. Most likely the two images are too dissimiliar.\n'.format(e.message))
        pass      
    if os.path.exists(ElastixResultFile) == True:    
        writer(out,'Elastix Registration Successfully Completed\n')
    #check to see if it was MHD instead
    elif os.path.exists(os.path.join(out, 'result.{}.mhd'.format((len(parameters)-1)))) == True:    
        ElastixResultFile = os.path.join(out, 'result.{}.mhd'.format((len(parameters)-1)))
        writer(out,'Elastix Registration Successfully Completed\n')
    else:
        writer(out, '\n***ERROR***Cannot find elastix result file\n: {}'.format(ElastixResultFile))
        return
        
    return ElastixResultFile, TransformParameterFile

#%%
def registration_fix(src, atlas, ann):

    #set paths        
    dst = os.path.join(src, 'clearmap_cluster_output') 
    auto = os.path.join(dst, 'autofluo_resampled.tif')     
    parameters = [os.path.join(src, 'clearmap_cluster/parameterfolder/Order1_Par0000affine.txt'),
                  os.path.join(src, 'clearmap_cluster/parameterfolder/Order2_Par0000bspline.txt')]
    
    #run transform
    masked_atlas_pth, masked_ann_pth = mask_similarity_transformed_atlas(auto, ann, atlas, dst, verbose=True)

    #elastix
    reg_dir = os.path.join(dst, 'elastix_auto_to_sim_atlas'); makedir(reg_dir)
    reg_dst = elastix_command_line_call(masked_atlas_pth, auto, reg_dir, parameters)
    return reg_dst
    