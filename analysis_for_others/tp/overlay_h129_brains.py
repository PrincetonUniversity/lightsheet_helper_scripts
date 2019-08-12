#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 09:37:20 2019

@author: wanglab
"""


import os, sys, numpy as np, pandas as pd, cv2, multiprocessing as mp
from skimage.external import tifffile
from skimage.morphology import ball
sys.path.append("/jukebox/wang/zahra/lightsheet_copy")
from tools.utils.io import makedir, load_dictionary, listdirfull, listall
from tools.registration.transform_list_of_points import modify_transform_files
from tools.registration.register import transformix_command_line_call
from tools.imageprocessing.orientation import fix_contour_orientation, fix_dimension_orientation
from tools.registration.transform_cell_counts import generate_transformed_cellcount, get_fullsizedims_from_kwargs, points_resample

def overlay_qc(args):  
    #unpacking this way for multiprocessing
    fld, folder_suffix, output_folder, verbose, doubletransform, make_volumes = args
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
            kwargs = load_dictionary(lightsheet_parameter_dictionary)
            vol = [xx for xx in kwargs["volumes"] if xx.ch_type == "cellch"][0].resampled_for_elastix_vol
            transformed_vol = os.path.join(dst, "transformed_volume"); makedir(transformed_vol)
            if not doubletransform:
                transformfiles = [os.path.join(fld, "elastix/TransformParameters.0.txt"), os.path.join(fld, 
                                  "elastix/TransformParameters.1.txt")]
                transformfiles = modify_transform_files(transformfiles, transformed_vol) #copy over elastix files
                transformix_command_line_call(vol, transformed_vol, transformfiles[-1])
            else:
                v=[xx for xx in kwargs["volumes"] if xx.ch_type == "cellch"][0]
                #sig to reg
                tps = [listall(os.path.dirname(v.ch_to_reg_to_atlas), "/TransformParameters.0")[0], 
                       listall(os.path.dirname(v.ch_to_reg_to_atlas), "/TransformParameters.1")[0]]
                #reg to atlas
                transformfiles = tps+[os.path.join(fld, "elastix/TransformParameters.0.txt"), 
                                      os.path.join(fld, "elastix/TransformParameters.1.txt")]
                transformfiles = modify_transform_files(transformfiles, transformed_vol) #copy over elastix files
                transformix_command_line_call(vol, transformed_vol, transformfiles[-1])
            

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
            r = 1
            selem = ball(r)[int(r/2)]
            cell_cnn = cell_cnn.astype("uint8")
            cell_cnn = np.asarray([cv2.dilate(cell_cnn[i], selem, iterations = 1) for i in range(cell_cnn.shape[0])])
            
            tarr=np.asarray(tarr)
            if len(badlist)>0: 
                print("{} errors in mapping with cell_cnn shape {}, each max dim {}, \npossibly due to a registration overshoot \
                      or not using double transform\n\n{}".format(len(badlist), cell_cnn.shape, np.max(tarr,0), badlist))
            merged = np.stack([cell_cnn, cell_reg, np.zeros_like(cell_reg)], -1)
            tifffile.imsave(os.path.join(transformed_vol, "merged.tif"), merged)#, compress=1)
            #out = np.concatenate([cell_cnn, cell_reg, ], 0)
        
        #####check at the resampled for elastix phase before transform...this mapping looks good...
        if make_volumes:
            #make zyx numpy arry
            zyx = dataframe[["z","y","x"]].values
            kwargs = load_dictionary(lightsheet_parameter_dictionary)
            vol = [xx for xx in kwargs["volumes"] if xx.ch_type =="cellch"][0]
            fullsizedimensions = get_fullsizedims_from_kwargs(kwargs) #don"t get from kwargs["volumes"][0].fullsizedimensions it"s bad! use this instead
            zyx = fix_contour_orientation(zyx, verbose=verbose, **kwargs) #now in orientation of resample
            zyx = points_resample(zyx, original_dims = fix_dimension_orientation(fullsizedimensions, **kwargs), 
                                  resample_dims = tifffile.imread(vol.resampled_for_elastix_vol).shape, verbose = verbose)[:, :3]
            
            #cell channel
            cell_ch = tifffile.imread(vol.resampled_for_elastix_vol)
            cell_cnn = np.zeros_like(cell_ch)
            tarr = []; badlist=[]
            for _zyx in zyx:
                z,y,x = [int(xx) for xx in _zyx]
                tarr.append([z,y,x])
                try:
                    cell_cnn[z,y,x] = 100
                except:
                    badlist.append([z,y,x])
            tarr=np.asarray(tarr)        
            merged = np.stack([cell_cnn, cell_ch, np.zeros_like(cell_ch)], -1)
            tifffile.imsave(os.path.join(transformed_vol, "resampled_merged.tif"), merged)#, compress=1)
            
    except Exception as e:
        print(e)
        with open(error_file, "a") as err_fl:
            err_fl.write("\n\n{} {}\n\n".format(fld, e))
                
#%%
if __name__ == "__main__":
    #goal is to transform cooridnates, voxelize based on number of cells and overlay with reigstered cell signal channel...
    output_folder = "/home/wanglab/Desktop/thal_transformed_points"; makedir(output_folder)
    qc_folder =  "/home/wanglab/Desktop/thal_qc"; makedir(qc_folder)
    error_file = "/home/wanglab/Desktop/thal_qc/errors.txt"
    
    #inputs
    input_folder = "/jukebox/wang/pisano/tracing_output/antero_4x/"
    folder_suffix = "3dunet_output/pooled_cell_measures"
#    brains = ["20180409_jg46_bl6_lob6a_04",
#             "20180608_jg75",
#             "20170204_tp_bl6_cri_1750r_03",
#             "20180608_jg72",
#             "20180416_jg56_bl6_lob8_04",
#             "20170116_tp_bl6_lob45_ml_11",
#             "20180417_jg60_bl6_cri_04",
#             "20180410_jg52_bl6_lob7_05",
#             "20170116_tp_bl6_lob7_1000r_10",
#             "20180409_jg44_bl6_lob6a_02",
#             "20180410_jg49_bl6_lob45_02",
#             "20180410_jg48_bl6_lob6a_01",
#             "20180612_jg80",
#             "20180608_jg71",
#             "20170212_tp_bl6_crii_1000r_02",
#             "20170115_tp_bl6_lob6a_rpv_03",
#             "20170212_tp_bl6_crii_2000r_03",
#             "20180417_jg58_bl6_sim_02",
#             "20170130_tp_bl6_sim_1750r_03",
#             "20170115_tp_bl6_lob6b_ml_04",
#             "20180410_jg50_bl6_lob6b_03",
#             "20170115_tp_bl6_lob6a_1000r_02",
#             "20170116_tp_bl6_lob45_500r_12",
#             "20180612_jg77",
#             "20180612_jg76",
#             "20180416_jg55_bl6_lob8_03",
#             "20170115_tp_bl6_lob6a_500r_01",
#             "20170130_tp_bl6_sim_rpv_01",
#             "20170204_tp_bl6_cri_1000r_02",
#             "20170212_tp_bl6_crii_250r_01",
#             "20180417_jg61_bl6_crii_05",
#             "20170116_tp_bl6_lob7_ml_08",
#             "20180409_jg47_bl6_lob6a_05"]
#    
    brains = ["20170410_tp_bl6_lob6a_ml_repro_01",
         "20160823_tp_bl6_cri_500r_02",
         "20180417_jg59_bl6_cri_03",
         "20170207_db_bl6_crii_1300r_02",
         "20160622_db_bl6_unk_01",
         "20161205_tp_bl6_sim_750r_03",
         "20180410_jg51_bl6_lob6b_04",
         "20170419_db_bl6_cri_rpv_53hr",
         "20170116_tp_bl6_lob6b_lpv_07",
         "20170411_db_bl6_crii_mid_53hr",
         "20160822_tp_bl6_crii_1500r_06",
         "20160920_tp_bl6_lob7_500r_03",
         "20170207_db_bl6_crii_rpv_01",
         "20161205_tp_bl6_sim_250r_02",
         "20161207_db_bl6_lob6a_500r_53hr",
         "20170130_tp_bl6_sim_rlat_05",
         "20170115_tp_bl6_lob6b_500r_05",
         "20170419_db_bl6_cri_mid_53hr",
         "20161207_db_bl6_lob6a_850r_53hr",
         "20160622_db_bl6_crii_52hr_01",
         "20161207_db_bl6_lob6a_50rml_53d5hr",
         "20161205_tp_bl6_lob45_1000r_01",
         "20160801_db_l7_cri_01_mid_64hr"]
    
    brains = [os.path.join(input_folder, xx) for xx in brains]
    input_list = [xx for xx in brains if os.path.exists(os.path.join(xx, folder_suffix))]
    verbose, doubletransform, make_volumes = True, True, True
    iterlst = [(fld, folder_suffix, output_folder, verbose, doubletransform, make_volumes) for fld in input_list]
    p = mp.Pool(6)
    p.map(overlay_qc, iterlst)
    
