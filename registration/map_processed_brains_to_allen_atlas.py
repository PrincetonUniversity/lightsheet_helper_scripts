#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:11:31 2019

@author: wanglab
"""

import os, tifffile, time, shutil, sys
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
import cv2
import numpy as np, subprocess as sp
from tools.imageprocessing.preprocessing import resample

#helpers
def makedir(dr):
    if not os.path.exists(os.path.dirname(dr)): os.mkdir(os.path.dirname(dr))
    if not os.path.exists(dr): os.mkdir(dr)
    
def listdirfull(pth):
    lst = [os.path.join(pth, xx) for xx in os.listdir(pth)]
    
    return lst

def resize_save_helper(outdr, resizefactor, brainname, zpln, ch, im):
    svloc = os.path.join(outdr, brainname+"_resized_ch"+ch)
    makedir(svloc)
    if len(im.shape) == 2:  #grayscale images
        y,x = im.shape
        xfct = round(x/resizefactor)
        yfct = round(y/resizefactor)
    im1 = cv2.resize(im, (xfct, yfct), interpolation=cv2.INTER_LINEAR) #interpolation=cv2.INTER_AREA) #
    tifffile.imsave(os.path.join(svloc, brainname + "_C"+ch+"_Z"+zpln+".tif"), im1.astype("uint16"))
    return svloc

def process_planes_from_fullsizedatafolder(volumes, job, cores, outdr, verbose=True):
    """Easy way to pull outfiles from fullsizedatafolder"""
    ############################inputs
    zpln = str(job).zfill(4)
    bitdepth = "uint16" #default to uint16 unless specified
    ####################################

    for vol in volumes:
        fl = [xx for xx in listdirfull(vol) if "Z{}".format(zpln) in xx]
        if len(fl) == 1: resize_save_helper(os.path.join(outdr, os.path.basename(os.path.dirname(os.path.dirname(vol)))), 
                                            3, os.path.basename(vol), zpln, vol[-2:], 
                                            tifffile.imread(fl[0]).astype(bitdepth))
    return

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

def change_transform_parameter_initial_transform(fl, initialtrnsformpth):
    '''
    (InitialTransformParametersFileName "NoInitialTransform")
    initialtrnsformpth = 'NoInitialTransform' or 'pth/to/transform.0.txt'
    '''
    fl1 = fl[:-5]+'0000.txt'
    
    with open(fl, 'r') as f, open(fl1, 'w') as new:
            for line in f:
                new.write('(InitialTransformParametersFileName "{}")\n'.format(initialtrnsformpth)) if 'InitialTransformParametersFileName' in line else new.write(line)
    
    #overwrite original transform file
    shutil.move(fl1, fl)
    return

def apply_transformix_and_register(reg_ch_resampledforelastix, sig_ch_resampledforelastix, sig_to_reg_out, reg_to_atl_out,
                                   AtlasFile, parameters, transformfile, resampled_zyx_dims):
  
    
    #run elastix on sig/inj channel -> reg channel (but don't register reg to itself)
    ElastixResultFile, TransformParameterFile = elastix_command_line_call(reg_ch_resampledforelastix, sig_ch_resampledforelastix, sig_to_reg_out, parameters)
    
    #copy transform paramters to set up transform series:
    [shutil.copy(os.path.join(reg_to_atl_out, xx), os.path.join(sig_to_reg_out, 'regtoatlas_'+xx)) for xx in os.listdir(reg_to_atl_out) if 'TransformParameters' in xx]
    
    #connect transforms by setting regtoatlas TP0's initial transform to sig->reg transform
    #might need to go backwards...
    reg_to_atlas_tps = [os.path.join(sig_to_reg_out, xx) for xx in os.listdir(sig_to_reg_out) if 'TransformParameters' in xx and 'regtoatlas' in xx]; reg_to_atlas_tps.sort() 
    sig_to_reg_tps = [os.path.join(sig_to_reg_out, xx) for xx in os.listdir(sig_to_reg_out) if 'TransformParameters' in xx and 'regtoatlas' not in xx]; sig_to_reg_tps.sort()

    #account for moving the reg_to_atlas_tps:
    [change_transform_parameter_initial_transform(reg_to_atlas_tps[xx+1], reg_to_atlas_tps[xx]) for xx in range(len(reg_to_atlas_tps)-1)]

    #now make the initialtransform of the first(0) sig_to_reg be the last's reg_to_atlas transform
    change_transform_parameter_initial_transform(reg_to_atlas_tps[0], sig_to_reg_tps[-1])
            
    #run transformix        
    sp.call(['transformix', '-in', sig_ch_resampledforelastix, '-out', sig_to_reg_out, '-tp', reg_to_atlas_tps[-1]])
    
    print(sig_to_reg_out,'\n   Transformix File Generated: {}'.format(sig_to_reg_out))
        
#%%
if __name__ == "__main__":
    
    #setup
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    outdr = "/jukebox/wang/zahra/h129_contra_vs_ipsi/reg_to_allen"
    src = "/jukebox/wang/pisano/tracing_output/antero_4x"
    
    animals = ['20180409_jg46_bl6_lob6a_04',
                 '20180608_jg75',
                 '20170204_tp_bl6_cri_1750r_03',
                 '20180608_jg72',
                 '20180416_jg56_bl6_lob8_04',
                 '20170116_tp_bl6_lob45_ml_11',
                 '20180417_jg60_bl6_cri_04',
                 '20180410_jg52_bl6_lob7_05',
                 '20170116_tp_bl6_lob7_1000r_10',
                 '20180409_jg44_bl6_lob6a_02',
                 '20180410_jg49_bl6_lob45_02',
                 '20180410_jg48_bl6_lob6a_01',
                 '20180612_jg80',
                 '20180608_jg71',
                 '20170212_tp_bl6_crii_1000r_02',
                 '20170115_tp_bl6_lob6a_rpv_03',
                 '20170212_tp_bl6_crii_2000r_03',
                 '20180417_jg58_bl6_sim_02',
                 '20170130_tp_bl6_sim_1750r_03',
                 '20170115_tp_bl6_lob6b_ml_04',
                 '20180410_jg50_bl6_lob6b_03',
                 '20170115_tp_bl6_lob6a_1000r_02',
                 '20170116_tp_bl6_lob45_500r_12',
                 '20180612_jg77',
                 '20180612_jg76',
                 '20180416_jg55_bl6_lob8_03',
                 '20170115_tp_bl6_lob6a_500r_01',
                 '20170130_tp_bl6_sim_rpv_01',
                 '20170204_tp_bl6_cri_1000r_02',
                 '20170212_tp_bl6_crii_250r_01',
                 '20180417_jg61_bl6_crii_05',
                 '20170116_tp_bl6_lob7_ml_08',
                 '20180409_jg47_bl6_lob6a_05']
    
    animal = animals[jobid]
    
    print(animal)
    
    dst = os.path.join(outdr, animal)
    
    volumes = [os.path.join((os.path.join(src, animal)), xx) for xx in 
               os.listdir(os.path.join(src, animal)) if "resampledforelastix.tif" in xx]; volumes.sort()

    #step3
    #reg to atlas
    
    print("\n\n registration channel to atlas")
    fx = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
    out_reg = os.path.join(os.path.join(outdr, animal), "reg_to_atl"); makedir(out_reg)
    
    mv = [xx for xx in volumes if "488" in xx and "ch00" in xx][0]
        
    params = ["/jukebox/wang/zahra/python/lightsheet_py3/parameterfolder/Order1_Par0000affine.txt", 
              "/jukebox/wang/zahra/python/lightsheet_py3/parameterfolder/Order2_Par0000bspline.txt"]
    
#    e_out_file, transformfile = elastix_command_line_call(fx, mv, out_reg, params, fx_mask=False, verbose=False)
    transformfile = os.path.join(out_reg, "TransformParameters.1.txt")
    secondary_registration = True
    resampled_zyx_dims = False
    
    #inj to reg
    inj_vol = [xx for xx in volumes if "488" in xx and "ch01" in xx][0]   
     
    out_inj = os.path.join(os.path.join(outdr, animal), "inj_to_reg"); makedir(out_inj)
    #appy transform
#    apply_transformix_and_register(mv, inj_vol, out_inj, out_reg,
#                                   fx, params, transformfile, resampled_zyx_dims)
    
    #cell to reg
    cell_vol = [xx for xx in volumes if "647" in xx and "ch00" in xx][0]   
     
    out_cell = os.path.join(os.path.join(outdr, animal), "cell_to_reg"); makedir(out_cell)
    #appy transform
    
    apply_transformix_and_register(mv, cell_vol, out_cell, out_reg,
                                   fx, params, transformfile, resampled_zyx_dims)
    
    
    ####### check to see if script finished due to an error
    if os.path.exists(out_reg)==False:
        print("****ERROR****GOTTEN TO END OF SCRIPT,\nTHIS ELASTIX OUTPUT FILE DOES NOT EXIST: {0} \n".format(out_reg))
#
#    #inj to reg
#    print("\n\n injection channel to registration channel")
#    out_inj = os.path.join(os.path.join(outdr, animal), "inj_to_reg"); makedir(out_inj)
#    fx = os.path.join(out_reg, "result.1.tif")
#    mv = [xx for xx in volumes if "488" in xx and "ch01" in xx][0]
#    
#    elastix_command_line_call(fx, mv, out_inj, params, fx_mask=False, verbose=False)
#    
#    #cell to reg
#    print("\n\n cell channel to registration channel")
#    out_cell = os.path.join(os.path.join(outdr, animal), "cell_to_reg"); makedir(out_cell)
#    fx = os.path.join(out_reg, "result.1.tif")
#    mv = [xx for xx in volumes if "647" in xx and "ch00" in xx][0]
#    
#    elastix_command_line_call(fx, mv, out_cell, params, fx_mask=False, verbose=False)