#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 16:51:59 2019

@author: wanglab
"""

import os, sys, numpy as np, pandas as pd, cv2, multiprocessing as mp
from skimage.external import tifffile
from skimage.morphology import ball
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.utils.io import makedir, load_kwargs, listdirfull, listall
from tools.registration.transform_list_of_points import modify_transform_files, point_transformix, unpack_pnts, create_text_file_for_elastix
from tools.registration.register import transformix_command_line_call
from tools.imageprocessing.orientation import fix_contour_orientation, fix_dimension_orientation
from tools.registration.transform_cell_counts import change_transform_parameter_initial_transform, points_resample

def get_fullsizedimensions(dct):
    """ works around param dict in case paths were missaved """    
    try:
        kwargs = load_kwargs(dct)
        vol = [xx for xx in kwargs["volumes"] if xx.ch_type =="cellch"][0]
        zf = len(listdirfull(vol.full_sizedatafld_vol, ".tif"))
        yf,xf = tifffile.imread(listdirfull(vol.full_sizedatafld_vol, "tif")[0]).shape
        fullsizedimensions = tuple((zf, yf, xf))
    except: #if param dict is messed up
        fsz = os.path.join(os.path.dirname(dct), "full_sizedatafld")
        vols = os.listdir(fsz); vols.sort()
        src = os.path.join(fsz, vols[len(vols)-1]) #hack - try to load param_dict instead?
        if not os.path.isdir(src): src = os.path.join(fsz, vols[len(vols)-2])     
        zf = len(listdirfull(src, ".tif"))
        yf,xf = tifffile.imread(listdirfull(src, "tif")[0]).shape
        fullsizedimensions = tuple((zf, yf, xf))
    
    return fullsizedimensions
    
def get_resampledvol_n_dimensions(dct):
    """ works around param dict in case paths were missaved """
    try:
        kwargs = load_kwargs(dct)
        vol = [xx for xx in kwargs["volumes"] if xx.ch_type =="cellch"][0]
        resampled_vol = vol.resampled_for_elastix_vol
        resampled_dims = tifffile.imread(resampled_vol).shape        
    except FileNotFoundError:
        fls = listdirfull(os.path.dirname(dct), ".tif"); fls.sort()
        resampled_vol = fls[-1] #will be the last one, bc of the 647 channel
        resampled_dims = tifffile.imread(resampled_vol).shape
        
    return resampled_dims, resampled_vol

def generate_transformed_cellcount(dataframe, dst, transformfiles, lightsheet_parameter_dictionary, verbose=False):
    """Function to take a csv file and generate an input to transformix
    
    Inputs
    ----------------
    dataframe = preloaded pandas dataframe
    dst = destination to save files
    transformfiles = list of all elastix transform files used, and in order of the original transform****
    lightsheet_parameter_file = .p file generated from lightsheet package
    """
    #set up locations
    transformed_dst = os.path.join(dst, "transformed_points"); makedir(transformed_dst)
    
    #make zyx numpy arry
    zyx = dataframe[["z","y","x"]].values
    
    #adjust for reorientation THEN rescaling, remember full size data needs dimension change releative to resample
    fullsizedimensions = get_fullsizedimensions(lightsheet_parameter_dictionary)
    kwargs = load_kwargs(lightsheet_parameter_dictionary)
     
    zyx = fix_contour_orientation(zyx, verbose=verbose, **kwargs) #now in orientation of resample
    resampled_dims, resampled_vol = get_resampledvol_n_dimensions(lightsheet_parameter_dictionary)
    
    zyx = points_resample(zyx, original_dims = fix_dimension_orientation(fullsizedimensions, 
            **kwargs), resample_dims = resampled_dims, verbose = verbose)[:, :3]
         
    #make into transformix-friendly text file
    pretransform_text_file = create_text_file_for_elastix(zyx, transformed_dst)
        
    #copy over elastix files
    transformfiles = modify_transform_files(transformfiles, transformed_dst) 
    change_transform_parameter_initial_transform(transformfiles[0], "NoInitialTransform")
   
    #run transformix on points
    points_file = point_transformix(pretransform_text_file, transformfiles[-1], transformed_dst)
    
    #convert registered points into structure counts
    converted_points = unpack_pnts(points_file, transformed_dst)   
    
    return converted_points
    
def overlay_qc(args):  
    #unpacking this way for multiprocessing
    fld, folder_suffix, output_folder, verbose, doubletransform, make_volumes = args
    
    #init error file
    error_file = os.path.join(os.path.join(output_folder, os.path.basename(fld)), "errors.txt")
    
    try:
        #get 3dunet cell dataframe csv file
        input_csv = listdirfull(os.path.join(fld, folder_suffix), ".csv")
        assert len(input_csv) == 1, "multiple csv files"
        dataframe = pd.read_csv(input_csv[0])
        
        #location to save out
        dst = os.path.join(output_folder, os.path.basename(fld)); makedir(dst)
    
        #EXAMPLE USING LIGHTSHEET - assumes marking centers in the "raw" full sized cell channel. This will transform those 
        #centers into "atlas" space (in this case the moving image)
        #in this case the "inverse transform has the atlas as the moving image in the first step, 
        #and the autofluorescence channel as the moving image in the second step 
        #NOTE - it seems that the registration of cell to auto is failing on occasion....thus get new files...
        ################################
        cell_inverse_folder = listdirfull(os.path.join(fld, "elastix_inverse_transform"), "cellch")[0]
        a2r = listall(cell_inverse_folder, "atlas2reg_TransformParameters"); a2r.sort()
        r2s = listall(cell_inverse_folder, "reg2sig_TransformParameters"); r2s.sort() #possibly remove

        #IMPORTANT. the idea is to apply cfos->auto->atlas
        transformfiles = r2s + a2r if doubletransform else a2r #might get rid of r2s
    
        lightsheet_parameter_dictionary = os.path.join(fld, "param_dict.p")
            
        converted_points = generate_transformed_cellcount(dataframe, dst, transformfiles, 
                                                          lightsheet_parameter_dictionary, verbose=verbose)
    
        #load and convert to single voxel loc
        zyx = np.asarray([str((int(xx[0]), int(xx[1]), int(xx[2]))) for xx in np.nan_to_num(np.load(converted_points))])
        from collections import Counter
        zyx_cnt = Counter(zyx)
        
        #check...
        if make_volumes:
            #manually call transformix
            kwargs = load_kwargs(lightsheet_parameter_dictionary)
            resampled_dims, resampled_vol = get_resampledvol_n_dimensions(lightsheet_parameter_dictionary)
            transformed_vol = os.path.join(dst, "transformed_volume"); makedir(transformed_vol)
            if not doubletransform:
                transformfiles = [os.path.join(fld, "elastix/TransformParameters.0.txt"), os.path.join(fld, 
                                  "elastix/TransformParameters.1.txt")]
                transformfiles = modify_transform_files(transformfiles, transformed_vol) #copy over elastix files
                transformix_command_line_call(resampled_vol, transformed_vol, transformfiles[-1])
            else:
                v=[xx for xx in kwargs["volumes"] if xx.ch_type == "cellch"][0]
                #sig to reg
                tps = [listall(os.path.dirname(v.ch_to_reg_to_atlas), "/TransformParameters.0")[0], 
                       listall(os.path.dirname(v.ch_to_reg_to_atlas), "/TransformParameters.1")[0]]
                #reg to atlas
                transformfiles = tps+[os.path.join(fld, "elastix/TransformParameters.0.txt"), 
                                      os.path.join(fld, "elastix/TransformParameters.1.txt")]
                transformfiles = modify_transform_files(transformfiles, transformed_vol) #copy over elastix files
                transformix_command_line_call(resampled_vol, transformed_vol, transformfiles[-1])
            

            #cell_registered channel
            cell_reg = tifffile.imread(os.path.join(transformed_vol, "result.tif"))
            tifffile.imsave(os.path.join(transformed_vol, "result.tif"), cell_reg, compress=1)
            cell_cnn = np.zeros_like(cell_reg)
            tarr = []; badlist=[]
            for zyx,v in zyx_cnt.items():
                z,y,x = [int(xx) for xx in zyx.replace("(","",).replace(")","").split(",")]
                tarr.append([z,y,x])
                try:
                    cell_cnn[z,y,x] = v*100
                except:
                    badlist.append([z,y,x])
                    
            #apply x y dilation
            r = 2
            selem = ball(r)[int(r/2)]
            cell_cnn = cell_cnn.astype("uint8")
            cell_cnn = np.asarray([cv2.dilate(cell_cnn[i], selem, iterations = 1) for i in range(cell_cnn.shape[0])])
            
            tarr=np.asarray(tarr)
            if len(badlist)>0: 
                print("{} errors in mapping with cell_cnn shape {}, each max dim {}, \npossibly due to a registration overshoot \
                      or not using double transform\n\n{}".format(len(badlist), cell_cnn.shape, np.max(tarr,0), badlist))
            merged = np.stack([cell_reg, cell_cnn, np.zeros_like(cell_reg)], -1)
            tifffile.imsave(os.path.join(transformed_vol, "merged.tif"), merged)#, compress=1)
            #out = np.concatenate([cell_cnn, cell_reg, ], 0)
        
            #####check at the resampled for elastix phase before transform
            #make zyx numpy arry
            zyx = dataframe[["z","y","x"]].values
            fullsizedimensions = get_fullsizedimensions(lightsheet_parameter_dictionary) 
            zyx = fix_contour_orientation(zyx, verbose=verbose, **kwargs) #now in orientation of resample
            zyx = points_resample(zyx, original_dims = fix_dimension_orientation(fullsizedimensions, **kwargs), 
                                  resample_dims = resampled_dims, verbose = verbose)[:, :3]
            
            #cell channel
            cell_ch = tifffile.imread(resampled_vol)
            cell_cnn = np.zeros_like(cell_ch)
            tarr = []; badlist=[]
            for _zyx in zyx:
                z,y,x = [int(xx) for xx in _zyx]
                tarr.append([z,y,x])
                try:
                    cell_cnn[z,y,x] = 100
                except:
                    badlist.append([z,y,x])
            tarr = np.asarray(tarr)        
            merged = np.stack([cell_ch, cell_cnn, np.zeros_like(cell_ch)], -1)
            tifffile.imsave(os.path.join(transformed_vol, "resampled_merged.tif"), merged)#, compress=1)
            
    except Exception as e:
        print(e)
        with open(error_file, "a") as err_fl:
            err_fl.write("\n\n{} {}\n\n".format(fld, e))
                

if __name__ == "__main__":
    
    #goal is to transform cooridnates, voxelize based on number of cells and overlay with reigstered cell signal channel...
    #set paths and make folders
    dst = "/jukebox/wang/zahra/prv/20190930"
    input_folder = "/jukebox/wang/pisano/tracing_output/retro_4x/"
    folder_suffix = "3dunet_output/pooled_cell_measures"

    output_folder = os.path.join(dst, "prv_transformed_cells"); makedir(output_folder)
    
    #inputs
    brains = ["20180322_jg_bl6f_prv_28"]
    
    brains = [os.path.join(input_folder, xx) for xx in brains]
    input_list = [xx for xx in brains if os.path.exists(os.path.join(xx, folder_suffix))]
    verbose, doubletransform, make_volumes = True, False, True
    iterlst = [(fld, folder_suffix, output_folder, verbose, doubletransform, make_volumes) for fld in input_list]
    p = mp.Pool(6)
    p.map(overlay_qc, iterlst)