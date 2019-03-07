#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:14:44 2018

@author: wanglab
"""

import os
from ClearMap.cluster.preprocessing import makedir, pth_update
from ClearMap.cluster.par_tools import load_kwargs
from ClearMap.parameter_file import set_parameters_for_clearmap
import ClearMap.IO as io
from ClearMap.Analysis.Statistics import thresholdPoints
from ClearMap.Alignment.Resampling import resamplePoints, resamplePointsInverse
from ClearMap.Alignment.Elastix import transformPoints
from ClearMap.Analysis.Voxelization import voxelize
from ClearMap.Analysis.Label import countPointsInRegions, labelToName
os.chdir('/jukebox/wang/zahra/lightsheet_copy')
from tools.registration.registration_using_similarity_mask import mask_similarity_transformed_atlas
from tools.utils.io import writer
import subprocess as sp
import numpy as np

ann = '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif'
atlas = '/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif'
src = '/jukebox/wang/Jess/lightsheet_output/201810_cfos/processed/dadult_pc_crus1_2'

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

def output_analysis_helper(threshold = (20, 900), row = (3,3), **params):
    '''
    Function to change elastix result directory before running 'step 6' i.e. point transformix to atlas.
    '''
    dct = pth_update(set_parameters_for_clearmap(**params))
    
    dct['RegistrationAlignmentParameter']["resultDirectory"] = os.path.join(params["outputdirectory"], 'clearmap_cluster_output/elastix_auto_to_sim_atlas')
    
    points, intensities = io.readPoints(dct['ImageProcessingParameter']["sink"]);
    
    #Thresholding: the threshold parameter is either intensity or size in voxel, depending on the chosen "row"
    #row = (0,0) : peak intensity from the raw data
    #row = (1,1) : peak intensity from the DoG filtered data
    #row = (2,2) : peak intensity from the background subtracted data
    #row = (3,3) : voxel size from the watershed
    points, intensities = thresholdPoints(points, intensities, threshold = threshold, row = row);
    #points, intensities = thresholdPoints(points, intensities, threshold = (20, 900), row = (2,2));
    io.writePoints(dct['FilteredCellsFile'], (points, intensities));
    
    # Transform point coordinates
    #############################
    points = io.readPoints(dct['CorrectionResamplingPointsParameter']["pointSource"]);
    points = resamplePoints(**dct['CorrectionResamplingPointsParameter']);
    points = transformPoints(points, transformDirectory = dct['CorrectionAlignmentParameter']["resultDirectory"], indices = False, resultDirectory = None);
    dct['CorrectionResamplingPointsInverseParameter']["pointSource"] = points;
    points = resamplePointsInverse(**dct['CorrectionResamplingPointsInverseParameter']);
    dct['RegistrationResamplingPointParameter']["pointSource"] = points;
    points = resamplePoints(**dct['RegistrationResamplingPointParameter']);
    points = transformPoints(points, transformDirectory = dct['RegistrationAlignmentParameter']["resultDirectory"], indices = False, resultDirectory = None);
    io.writePoints(dct['TransformedCellsFile'], points);
   
    # Heat map generation
    #####################
    points = io.readPoints(dct['TransformedCellsFile'])
    intensities = io.readPoints(dct['FilteredCellsFile'][1])
    
    #Without weigths:
    vox = voxelize(points, dct['AtlasFile'], **dct['voxelizeParameter']);
    if not isinstance(vox, basestring):
      io.writeData(os.path.join(dct['OutputDirectory'], 'cells_heatmap.tif'), vox.astype('int32'));
    
    #With weigths from the intensity file (here raw intensity):
    dct['voxelizeParameter']["weights"] = intensities[:,0].astype(float);
    vox = voxelize(points, dct['AtlasFile'], **dct['voxelizeParameter']);
    if not isinstance(vox, basestring):
      io.writeData(os.path.join(dct['OutputDirectory'], 'cells_heatmap_weighted.tif'), vox.astype('int32'));
   
    #Table generation:
    ##################
    #With integrated weigths from the intensity file (here raw intensity):
    ids, counts = countPointsInRegions(points, labeledImage = dct['AnnotationFile'], intensities = intensities, intensityRow = 0);
    table = np.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
    table["id"] = ids;
    table["counts"] = counts;
    table["name"] = labelToName(ids);
    io.writeTable(os.path.join(dct['OutputDirectory'], 'Annotated_counts_intensities.csv'), table);
    
    #Without weigths (pure cell number):
    ids, counts = countPointsInRegions(points, labeledImage = dct['AnnotationFile'], intensities = None);
    table = np.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
    table["id"] = ids;
    table["counts"] = counts;
    table["name"] = labelToName(ids);
    io.writeTable(os.path.join(dct['OutputDirectory'], 'Annotated_counts.csv'), table);
    
    print ('Analysis Completed')    
    
    return

#%%
def registration_fix(src, atlas, ann):

    #set paths        
    kwargs = load_kwargs(src)
    dst = os.path.join(src, 'clearmap_cluster_output') 
    auto = os.path.join(dst, 'autofluo_resampled.tif')     
    parameters = [os.path.join(src, 'clearmap_cluster/parameterfolder/Order1_Par0000affine.txt'),
                  os.path.join(src, 'clearmap_cluster/parameterfolder/Order2_Par0000bspline.txt')]
    
    #run transform
    masked_atlas_pth, masked_ann_pth = mask_similarity_transformed_atlas(auto, ann, atlas, dst, verbose=True)

    #elastix
    reg_dir = os.path.join(dst, 'elastix_auto_to_sim_atlas'); makedir(reg_dir)
    reg_dst = elastix_command_line_call(masked_atlas_pth, auto, reg_dir, parameters)
    
    #'step 6' - point transformix
    output_analysis_helper(**kwargs)
    
    return reg_dst
    