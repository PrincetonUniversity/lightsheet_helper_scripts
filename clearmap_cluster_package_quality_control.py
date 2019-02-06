#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 13:55:49 2019

@author: wanglab
"""

import os, shutil, numpy as np
from skimage.external import tifffile
import skimage
import matplotlib.gridspec as gridspec
from scipy.ndimage.interpolation import zoom
from collections import Counter
os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from tools.utils.io import makedir, load_dictionary, load_np, listdirfull, load_memmap_arr, listall
from tools.registration.transform_list_of_points import create_text_file_for_elastix, point_transformix, modify_transform_files, unpack_pnts
from tools.registration.register import change_transform_parameter_initial_transform
from tools.imageprocessing.orientation import fix_dimension_orientation
from tools.registration.transform_cell_counts import points_resample, get_fullsizedims_from_kwargs
os.chdir("/jukebox/wang/zahra/clearmap_cluster_copy")
from ClearMap.cluster.utils import load_kwargs


def generate_transformed_cellcount(dataframe, dst, transformfiles, dct, verbose=False):
    '''Function to take a csv file and generate an input to transformix
    
    Inputs
    ----------------
    dataframe = preloaded pandas dataframe
    dst = destination to save files
    transformfiles = list of all elastix transform files used, and in order of the original transform****
    lightsheet_parameter_file = .p file generated from lightsheet package
    '''
    #set up locations
    transformed_dst = os.path.join(dst, 'transformed_points'); makedir(transformed_dst)
    
    #make zyx numpy arry
    zyx = dataframe
    
    #adjust for reorientation THEN rescaling, remember full size data needs dimension change releative to resample
    kwargs = load_dictionary(dct)
    zyx = zyx[:, [2, 1, 0]] #fix orientation  
    fullsizedimensions = get_fullsizedims_from_kwargs(kwargs) #don't get from kwargs['volumes'][0].fullsizedimensions it's bad! use this instead
    zyx = points_resample(zyx, original_dims = (fullsizedimensions[2], fullsizedimensions[1], fullsizedimensions[0]), 
                          resample_dims = tifffile.imread(os.path.join(fld, "clearmap_cluster_output/cfos_resampled.tif")).shape, 
                          verbose = verbose)[:, :3]
   
    #make into transformix-friendly text file
    pretransform_text_file = create_text_file_for_elastix(zyx, transformed_dst)
        
    #copy over elastix files
    transformfiles = modify_transform_files(transformfiles, transformed_dst) 
    change_transform_parameter_initial_transform(transformfiles[0], 'NoInitialTransform')
   
    #run transformix on points
    points_file = point_transformix(pretransform_text_file, transformfiles[-1], transformed_dst)
    
    #convert registered points into structure counts
    converted_points = unpack_pnts(points_file, transformed_dst)   
    
    return converted_points
#%%
if __name__ == '__main__':

    #set up
    dst = '/jukebox/wang/Jess/lightsheet_output/201812_development/forebrain/qc'; makedir(dst)
    lst = [xx for xx in listdirfull('/jukebox/wang/Jess/lightsheet_output/201812_development/forebrain/processed')]
    transform = 'all';#both for regwatlas, and only affine for sig adn reg #'all', 'single': don't consider reg with sig at all
    verbose = True
    #fld = '/home/wanglab/wang/pisano/tracing_output/antero_4x/20170130_tp_bl6_sim_1750r_03'    
    #loop
    for fld in lst:
        try:
            print fld
            kwargs = load_kwargs(fld)        
            dst0 = os.path.join(dst, os.path.basename(fld)); makedir(dst0)
            dst1 = os.path.join(dst0, 'elastix'); makedir(dst1)
            
            #####check cell detection (modeled from lightsheet/tools/registration/transform_cell_counts)
            #3dunet cell dataframe
            dataframe = load_np(os.path.join(fld, 'clearmap_cluster_output/cells.npy'))
        
            #assumes marking centers in the 'raw' full sized cell channel. This will transform those centers into "atlas" space (in this case the moving image)
            #in this case the "inverse transform has the atlas as the moving image in the first step, and the autofluorescence channel as the moving image in the second step 
            r2s0 = os.path.join(fld, "clearmap_cluster_output/elastix_cfos_to_auto/TransformParameters.0.txt")
            r2s1 = os.path.join(fld, "clearmap_cluster_output/elastix_cfos_to_auto/TransformParameters.1.txt")
            a2r0 = os.path.join(fld, "clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters.0.txt")
            a2r1 = os.path.join(fld, "clearmap_cluster_output/elastix_auto_to_atlas/TransformParameters.1.txt")
            if transform == 'all':
                transformfiles = [r2s0, r2s1, a2r0, a2r1]
            elif transform == 'single':
                transformfiles = [a2r0, a2r1]
            elif transform == 'affine_only_reg_to_sig':    
                transformfiles = [r2s0, a2r0, a2r1]
            transformfiles = modify_transform_files(transformfiles, dst = dst1)
            
            #convert points
            converted_points = generate_transformed_cellcount(dataframe, dst1, transformfiles, 
                                                              dct = os.path.join(fld, 'param_dict.p'), verbose=verbose)
            
            #load and convert to single voxel loc
            zyx = np.asarray([str((int(xx[0]), int(xx[1]), int(xx[2]))) for xx in load_np(converted_points)])
            zyx_cnt = Counter(zyx)
                                
            #atlas
            atl = tifffile.imread("/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif")
            atl_cnn = np.zeros_like(atl)
            errors = []
            for zyx,v in zyx_cnt.iteritems():
                z,y,x = [int(xx) for xx in zyx.replace('(','',).replace(')','').split(',')]
                try:
                    atl[z,y,x] = v*100
                except Exception, e:
                    print e
                    errors.append(e)
            if len(errors)>0:
                with open(os.path.join(dst, '{}_errors.txt'.format(os.path.basename(fld))), 'a') as flll:
                    for err in errors:
                        flll.write(str(err)+'\n')
                    flll.close()
            merged = np.stack([atl_cnn, atl, np.zeros_like(atl)], -1)
            #reorient to horizontal
            merged = np.swapaxes(merged, 0,2)
            tifffile.imsave(os.path.join(dst, '{}_points_merged.tif'.format(os.path.basename(fld))), merged)
            shutil.rmtree(dst0)
            
            
            ####CHECK REGS       
            #get regs
            print fld, '\ngenerating registration metrics, loading....'        
            #get elastix paths
            ereg = os.path.join(regvol.elastixfld, 'result.1.tif')
            einj = injvol.ch_to_reg_to_atlas
            ecell = cellvol.ch_to_reg_to_atlas
            
            #get inverse
            a2rinj = [xx for xx in listall(injvol.inverse_elastixfld, '_atlas2reg/result.1.tif') if 'injch' in xx][0]
            r2sinj = [xx for xx in listall(injvol.inverse_elastixfld, '_reg2sig/result.1.tif') if 'injch' in xx][0]
            a2rcell = [xx for xx in listall(cellvol.inverse_elastixfld, '_atlas2reg/result.1.tif') if 'cellch' in xx][0]
            r2scell = [xx for xx in listall(cellvol.inverse_elastixfld, '_reg2sig/result.1.tif') if 'cellch' in xx][0]
            
            #load
            #vols_pths = [ereg, einj, ecell, a2rinj, r2sinj, a2rcell, r2scell]
            #names = ['auto elastix', 'inj elastix', 'cell elastix', 'atlas to auto', 'auto to inj', 'atlas to auto', 'auto to cell']
            tvols_pths = [ereg, einj, ecell, a2rinj, r2sinj, r2scell]
            tnames = ['auto\nelastix', 'inj\nelastix', 'cell\nelastix', 'atlas to\nauto', 'auto to\ninj', 'auto to\ncell']
            vols_pths=[]
            names=[]
            #this ensure things don't fail...
            for i in range(len(tvols_pths)):
                if os.path.exists(tvols_pths[i]):
                    vols_pths.append(tvols_pths[i])
                    names.append(tnames[i])
                else:
                    print('\n\nMissing {}'.format(tnames[i]))
            vols = [np.swapaxes(tifffile.imread(xx), 0, 2) for xx in vols_pths]
            
            #make figure
            fig = plt.figure(figsize=(6,8))
            rows = len(names)
            cols = 5
            tick = 0
            gs1 = gridspec.GridSpec(rows,cols)
            gs1.update(wspace=0.0000025, hspace=0.000005)
            adjust_exposure = True
            for row in range(rows):
                im = np.copy(vols[row])#*100
                step = im.shape[0] / cols
                for col in range(cols):
                    ax = plt.subplot(gs1[tick])
                    tick+=1; print tick
                    #ax = plt.subplot(rows,cols,tick)        
                    immax = np.max(im[col*step:(col+1)*step], 0)
                    if adjust_exposure: immax = skimage.exposure.equalize_hist(immax, nbins=20000)            
                    ax.imshow(immax, cmap='gray')
                    ax.set_yticklabels([])
                    ax.set_xticklabels([])
                    if col == 0: ax.set_ylabel(names[row]).set_rotation(45)
            fig.suptitle(os.path.basename(fld));                
            fig.subplots_adjust(top=0.88, wspace=0, hspace=0)
            plt.savefig(os.path.join(dst, os.path.basename(fld)+'_registration_qc.pdf'), dpi=300, transparent=True)
            plt.close()
            
            if False:
                #overlay of cells - non used, bitdepth issues
                #load
                #vols_pths = [ereg, einj, ecell, a2rinj, r2sinj, a2rcell, r2scell]
                #names = ['auto elastix', 'inj elastix', 'cell elastix', 'atlas to auto', 'auto to inj', 'atlas to auto', 'auto to cell']
                tvols_pths = [ereg, einj, ecell, a2rinj, r2sinj, r2scell]
                tnames = ['auto\nelastix', 'inj\nelastix', 'cell\nelastix', 'atlas to\nauto', 'auto to\ninj', 'auto to\ncell']
                vols_pths=[]
                names=[]
                #this ensure things don't fail...
                for i in range(len(tvols_pths)):
                    if os.path.exists(tvols_pths[i]):
                        vols_pths.append(tvols_pths[i])
                        names.append(tnames[i])
                    else:
                        print('\n\nMissing {}'.format(tnames[i]))
                vols = [np.swapaxes(tifffile.imread(xx), 0, 2) for xx in vols_pths]
                #done this way to prevent swapping of axes
                cnn_vol = tifffile.imread(os.path.join(dst, '{}_points_merged.tif'.format(os.path.basename(fld))))
                for i in range(3):
                    cnn_vol[i] = skimage.exposure.equalize_hist(cnn_vol[i], nbins=45000)
                vols.append(cnn_vol)
                names.append('CNN\noverlay')
                
                #make figure
                fig = plt.figure(figsize=(6,8))
                rows = len(names)
                cols = 5
                tick = 0
                gs1 = gridspec.GridSpec(rows,cols)
                gs1.update(wspace=0.0000025, hspace=0.000005)
                adjust_exposure = True
                for row in range(rows):
                    name = names[row]
                    im = np.copy(vols[row])#*100
                    step = im.shape[0] / cols
                    for col in range(cols):
                        ax = plt.subplot(gs1[tick])
                        tick+=1; print tick
                        #ax = plt.subplot(rows,cols,tick)
                        if name == 'CNN\noverlay':
                            immax = np.swapaxes(np.swapaxes(np.asarray([np.max(im[col*step:(col+1)*step,:,:,i], 0) for i in range(3)]),0,1),1,2)
                            if adjust_exposure: immax = np.swapaxes(np.swapaxes(np.asarray([skimage.exposure.equalize_hist(immax[:,:,i], nbins=20000) for i in range(3)]),0,1),1,2)
                        else:
                            immax = np.max(im[col*step:(col+1)*step], 0)
                            if adjust_exposure: immax = skimage.exposure.equalize_hist(immax, nbins=20000)            
                            
                        
                        ax.imshow(immax, cmap='gray')
                        ax.set_yticklabels([])
                        ax.set_xticklabels([])
                        if col == 0: ax.set_ylabel(name).set_rotation(45)
                fig.suptitle(os.path.basename(fld));                
                fig.subplots_adjust(top=0.88, wspace=0, hspace=0)
                plt.savefig(os.path.join(dst, os.path.basename(fld)+'_registration_qc.pdf'), dpi=300, transparent=True)
                plt.close()
        except Exception, e:
            print e
        
    