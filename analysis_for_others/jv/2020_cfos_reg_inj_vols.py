#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:54:26 2020

@author: wanglab
"""

import os, tifffile, time, shutil
import cv2
import numpy as np, subprocess as sp
from scipy.ndimage import zoom

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

#%%
if __name__ == "__main__":
    
    #setup
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    outdr = "/jukebox/wang/Jess/lightsheet_output/202002_cfos/injection/processed"
    src = "/jukebox/wang/Jess/lightsheet_output/202002_cfos/processed"
    
    animals = ["an10_crus1_lat", "an10_vecctrl_ymaze", "an11_crus1_lat",
       "an12_crus1_lat", "an13_crus1_lat", "an14_crus1_lat",
       "an15_crus1_lat", "an16_crus1_lat", "an17_crus1_lat",
       "an18_crus1_lat", "an19_crus1_lat", "an1_crus1_lat", "an1_saline",
       "an1_vecctrl_ymaze", "an20_crus1_lat", "an21_crus1_lat",
       "an22_crus1_lat", "an23_crus1_lat", "an24_crus1_lat",
       "an25_crus1_lat", "an26_crus1_lat", "an27_crus1_lat",
       "an28_crus1_lat", "an29_crus1_lat", "an2_crus1_lat", "an2_saline",
       "an2_vecctrl_ymaze", "an30_crus1_lat", "an31_crus1_lat",
       "an32_crus1_lat", "an33_crus1_lat", "an34_crus1_lat",
       "an3_crus1_lat", "an3_saline", "an3_vecctrl_ymaze",
       "an4_crus1_lat", "an4_saline", "an4_vecctrl_ymaze", "an5_cno",
       "an5_crus1_lat", "an5_vecctrl_ymaze", "an6_cno", "an6_crus1_lat",
       "an6_vecctrl_ymaze", "an7_cno", "an7_crus1_lat",
       "an7_vecctrl_ymaze", "an8_cno", "an8_crus1_lat",
       "an8_vecctrl_ymaze", "an9_crus1_lat", "an9_vecctrl_ymaze"]
    
    animal = animals[jobid]
    
    dst = os.path.join(outdr, animal)
    
    volumes = [os.path.join((os.path.join(src, animal)), "full_sizedatafld/"+xx) for xx in 
              os.listdir(os.path.join(os.path.join(src, animal), "full_sizedatafld")) if not xx[-3:] == "txt"]
    
    #part of step1    
    start = time.time()
    
    for job in range(900):
        process_planes_from_fullsizedatafolder(volumes, job, 12, outdr, verbose=True)
        
    print("\n\ntook {} minutes".format((time.time()-start)/60))    
    
    #step2
    
    volumes = listdirfull(dst); volumes.sort()
    
    for vol in volumes:
        plns = listdirfull(vol); plns.sort()
        y,x = tifffile.imread(plns[0]).shape
        #set destination
        memmap_dst = os.path.join(dst, os.path.join(os.path.basename(vol) + ".npy"))
        resz = np.lib.format.open_memmap(memmap_dst, dtype = 'uint16', mode = 'w+', shape = (len(plns), y, x))
        for i, pln in enumerate(plns):
            resz[i] = tifffile.imread(pln)
            if i%50 == 0: resz.flush()
        resz = np.transpose(resz, [2, 1, 0]) #sagittal
        img = os.path.join(dst, os.path.join(os.path.basename(vol)+".tif"))
        tifffile.imsave(img, resz)
        #delete unnecessary things once we have the image
        os.remove(memmap_dst)
        shutil.rmtree(vol)

    #step3
    #reg to atlas
    fx = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
    out = os.path.join(os.path.join(outdr, animal), "elastix"); makedir(out)
    
    mv = [os.path.join(dst, xx) for xx in os.listdir(dst) if "488" in xx][0]
    resz_shp = (702, 832, 457)
    img = tifffile.imread(mv)
    #resize
    resmpld = zoom(img, (resz_shp[0]/img.shape[0], resz_shp[1]/img.shape[1], resz_shp[2]/img.shape[2]), order = 3)
    
    #save out, overwrite
    tifffile.imsave(mv, resmpld)
#    
    params = ["/jukebox/wang/zahra/python/lightsheet_py3/parameterfolder/Order1_Par0000affine.txt", 
              "/jukebox/wang/zahra/python/lightsheet_py3/parameterfolder/Order2_Par0000bspline.txt"]
#    
    elastix_command_line_call(fx, mv, out, params, fx_mask=False, verbose=False)
#    
    #inj to reg
    out = os.path.join(os.path.join(outdr, animal), "elastix")
    fx = os.path.join(out, "result.1.tif")
    mv = [os.path.join(dst, xx) for xx in os.listdir(dst) if "647" in xx][0]
    resz_shp = (702, 832, 457)
    img = tifffile.imread(mv)
    #resize
    resmpld = zoom(img, (resz_shp[0]/img.shape[0], resz_shp[1]/img.shape[1], resz_shp[2]/img.shape[2]), order = 3)
    
    #save out, overwrite
    tifffile.imsave(mv, resmpld)
    
    out = os.path.join(out, os.path.basename(mv)[:-17]); makedir(out)
    elastix_command_line_call(fx, mv, out, params, fx_mask=False, verbose=False)