#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 17:18:52 2019

@author: wanglab
"""


import os, tifffile, time, shutil
import cv2
import numpy as np, subprocess as sp
from scipy.ndimage import zoom
from tools.utils.io import writer
from tools.imageprocessing.preprocessing import flatten_stitcher, stitcher, resize_save, saver

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

def process_planes(job, cores, compression, volumes, verbose=True):
    '''Function for Slurm Array job to process single zplane; this could be parallelized'''
    ############################inputs
    zpln = str(job).zfill(4)
    resizefactor = 3
    bitdepth = 'uint16' #default to uint16 unless specified
    ####################################
    
    for i, vol in enumerate(volumes):
        if i==0: writer(vol.full_sizedatafld, 'Processing zpln {}...'.format(zpln), flnm='process.txt', verbose=verbose)
        try:
            dct=vol.zdct[zpln] #dictionary of files for single z plane
        except KeyError:
            return 'ArrayJobID/SF exceeds number of planes'

        #handle raw vs nonraw
        if vol.raw: stitchdct = flatten_stitcher(cores, vol.outdr, vol.tiling_overlap, vol.xtile, vol.ytile, zpln, dct, vol.lightsheets, **kwargs)
        if not vol.raw: stitchdct = stitcher(cores, vol.outdr, vol.tiling_overlap, vol.xtile, vol.ytile, zpln, dct, **kwargs)

        #save for eventual compression
        saver(cores, stitchdct, vol.full_sizedatafld, vol.brainname, zpln, compression, bitdepth)

        #####################Cell detect ###CAN ADD A INJ CHANNEL HERE
        #resize
        resize_save(cores, stitchdct, vol.outdr, resizefactor, vol.brainname, zpln, bitdepth) ##can return location of cell channel but not interested in this; if so remove [0]
        writer(vol.full_sizedatafld, '\n   completed zpln {}, channel {}, channel_type {}, nm {}'.format(zpln, vol.channel, vol.ch_type, os.path.basename(vol.full_sizedatafld_vol)), flnm='process.txt', verbose=verbose)

    #log and return
    writer(vol.full_sizedatafld, '\n   ...completed all volumes of zpln {}\n'.format(zpln), flnm='process.txt', verbose=verbose)   
    
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

if __name__ == "__main__":
    
    #setup
    jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    
    outdr = "/jukebox/wang/pisano/figures/cfos_injection"
    src = "/jukebox/wang/pisano/tracing_output/cfos/201701_cfos/injection_site"
    
    animals = ['201701_mk06',
    '201701_tp02',
     '201701_tp01',
     '201701_tp05',
     '201701_tpbe',
     '201701_tp07',
    '201701_tpal',
     '201701_tp06',
     '201701_mk01',
     '201701_tp08',
     '201701_mk05',
     '201701_mk07',
     '201701_mk03',
     '201701_tp09',
     '201701_mk11',
     '201701_mk02',
     '201701_mk04',
     '201701_mk08',
     '201701_mk10']
    
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
    params = ["/jukebox/wang/zahra/lightsheet_copy/parameterfolder/Order1_Par0000affine.txt", 
              "/jukebox/wang/zahra/lightsheet_copy/parameterfolder/Order2_Par0000bspline.txt"]
#    
#    elastix_command_line_call(fx, mv, out, params, fx_mask=False, verbose=False)
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