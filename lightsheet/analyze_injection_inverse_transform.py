#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:07:15 2018

@author: wanglab
"""

import os, sys, numpy as np
import tifffile
import matplotlib as mpl
from tools.imageprocessing.orientation import fix_orientation
from tools.registration.transform import count_structure_lister, transformed_pnts_to_allen_helper_func
from tools.analysis.analyze_injection import orientation_crop_check, find_site
from tools.registration.register import make_inverse_transform, point_transform_due_to_resizing, point_transformix
from tools.utils.io import load_kwargs, listdirfull, makedir
from natsort import natsorted
from collections import Counter
import SimpleITK as sitk, pandas as pd
import matplotlib.pyplot as plt; plt.ion()
from scipy.ndimage import zoom

if __name__ == '__main__':
    
    #check if reorientation is necessary
    src = '/home/wanglab/LightSheetTransfer/tp/20200930_17_32_58_hsv_36hr_7/Ex_642_Em_2/downsized_for_atlas.tif'
    src = orientation_crop_check(src, axes = ('2','1','0'), crop = '[:,500:,:]') #'[:,390:,:]'
    
    #optimize detection parameters for inj det
    optimize_inj_detect(src, threshold =9, filter_kernel = (5,5,5))

    #run
    #suggestion: save_individual=True,
    inputlist = [
        '/home/wanglab/LightSheetTransfer/tp/20200930_17_32_58_hsv_36hr_7/Ex_642_Em_2/downsized_for_atlas.tif',
        '/home/wanglab/LightSheetTransfer/tp/20201001_10_57_49_hsv_36h_6/Ex_642_Em_2/downsized_for_atlas.tif',
        '/home/wanglab/LightSheetTransfer/tp/20201001_17_13_35_hsv_28h_2/Ex_642_Em_2/downsized_for_atlas.tif',
        '/home/wanglab/LightSheetTransfer/tp/20201001_15_39_26_hsv_28h_4/Ex_642_Em_2/downsized_for_atlas.tif',
        '/home/wanglab/LightSheetTransfer/tp/PRV_50hr-019/Ex_642_Em_2/downsized_for_atlas.tif',
        '/home/wanglab/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_PRV/natneuroreviews_tompisano_PRV_36hr-015/imaging_request_1/output/processing_request_1/resolution_4x/Ex_642_Em_2/downsized_for_atlas.tif',
        '/home/wanglab/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_PRV/natneuroreviews_tompisano_PRV_28hr-011/imaging_request_1/output/processing_request_1/resolution_4x/Ex_642_Em_2/downsized_for_atlas.tif',
        '/home/wanglab/wang/pisano/tracing_output/bl6_ts/20150804_tp_bl6_ts04/20150804_tp_bl6_ts04_555_z3um_70msec_3hfds_resized_ch00_resampledforelastix.tif',
        '/home/wanglab/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_CTB/natneuroreviews_tompisano_CTB-001/imaging_request_1/output/processing_request_1/resolution_4x/Ex_561_Em_1/downsized_for_atlas.tif',
        '/home/wanglab/LightSheetData/lightserv/jverpeut/natneuroreviews_tompisano_CTB/natneuroreviews_tompisano_CTB-002/imaging_request_1/output/processing_request_1/resolution_4x/Ex_561_Em_1/downsized_for_atlas.tif'
        ]
    
    kwargs = {'inputlist': inputlist,
          'channel': '01',
          'channel_type': 'injch',
          'filter_kernel': (5,5,5), #rbdg = (5,5,5)
          'threshold': 10, #rbdg = 10 NOTE: thresholding is different than analyze_injection.py
          'num_sites_to_keep': 1,
          'injectionscale': 45000, 
          'imagescale': 2,
          'reorientation': ('2','0','1'),
          'crop': '[:,500:,:]', #limits injection site search to cerebellum
          'dst': '/home/wanglab/Desktop/inj',
          'save_individual': True, 
          'colormap': 'plasma', 
          'atlas': '/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif',
          'annotation': '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif',
          'id_table': '/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx'#if using allen for registration, need to specify allen id table - path defaults to PMA id tables
        }
    
    df = pool_injections_inversetransform(**kwargs)
    
#%%
def pool_injections_inversetransform(**kwargs):
    '''Function to pool several injection sites. 
    Assumes that the basic registration AND inverse transform using elastix has been run. 
    If not, runs inverse transform. Additions to analyze_injection.py and pool_injections_for_analysis().

    Inputs
    -----------
    kwargs:
      'inputlist': inputlist, #list of folders generated previously from software
      'channel': '01', 
      'channel_type': 'injch',
      'filter_kernel': (5,5,5), #gaussian blur in pixels (if registered to ABA then 1px likely is 25um)
      'threshold': 10 (int, value to use for thresholding, this value represents the number of stand devs above the mean of the gblurred image)
      'num_sites_to_keep': #int, number of injection sites to keep, useful if multiple distinct sites
      'injectionscale': 45000, #use to increase intensity of injection site visualizations generated - DOES NOT AFFECT DATA
      'imagescale': 2, #use to increase intensity of background  site visualizations generated - DOES NOT AFFECT DATA
      'reorientation': ('2','0','1'), #use to change image orientation for visualization only
      'crop': #use to crop volume, values below assume horizontal imaging and sagittal atlas
                False
                cerebellum: '[:,600:,:]'
                caudal midbrain: '[:,300:415,:]'
                midbrain: '[:,215:415,:]'
                thalamus: '[:,215:345,:]'
                anterior cortex: '[:,:250,:]'
      
      'dst': '/home/wanglab/Downloads/test', #save location
      'save_individual': True, #optional to save individual images, useful to inspect brains, which you can then remove bad brains from list and rerun function
      'colormap': 'plasma', 
      'atlas': '/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif', #whole brain atlas
      
      Optional:
          ----------
          'save_array': path to folder to save out numpy array per brain of binarized detected site
          'save_tif': saves out tif volume per brain of binarized detected site
          'dpi': dots per square inch to save at

      Returns
      ----------------count_threshold
      a pooled image consisting of max IP of reorientations provide in kwargs.
      a list of structures (csv file) with pixel counts, pooling across brains.
      if save individual will save individual images, useful for inspection and/or visualization
    '''
    
    inputlist = kwargs['inputlist']
    dst = kwargs['dst']; makedir(dst)
    injscale = kwargs['injectionscale'] if 'injectionscale' in kwargs else 1
    imagescale = kwargs['imagescale'] if 'imagescale' in kwargs else 1
    axes = kwargs['reorientation'] if 'reorientation' in kwargs else ('0','1','2')
    cmap = kwargs['colormap'] if 'colormap' in kwargs else 'plasma'
    save_array = kwargs['save_array'] if 'save_array' in kwargs else False
    save_tif = kwargs['save_tif'] if 'save_tif' in kwargs else False
    num_sites_to_keep = kwargs['num_sites_to_keep'] if 'num_sites_to_keep' in kwargs else 1
    ann = sitk.GetArrayFromImage(sitk.ReadImage(kwargs['annotation']))
    #if kwargs['crop']: (from original analyze injection function, no functionality here if points file exist)
    #    ann = eval('ann{}'.format(kwargs['crop']))
    nonzeros = []
    #not needed as mapped points from point_transformix used
    id_table = kwargs['id_table'] if 'id_table' in kwargs else '/jukebox/wang/zahra/lightsheet_copy/supp_files/ls_id_table_w_voxelcounts.xlsx'
    id_table = pd.read_excel(id_table)
    
    for i in range(len(inputlist)): #to iteratre through brains
        outdr = dst
        
        im = tifffile.imread(inputlist[i]) #load inj_vol as numpy array
        if kwargs['crop']: im = eval('im{}'.format(kwargs['crop']))#; print im.shape
        
        print('  loading:\n     {}'.format(pth))
        #run find site function to segment inj site using non-registered resampled for elastix volume - pulled directly from tools.registration.register.py and tools.analysis.analyze_injection.py
        array = find_site(im, thresh=kwargs['threshold'], filter_kernel=kwargs['filter_kernel'], num_sites_to_keep = num_sites_to_keep)*injscale
        if save_array: np.save(os.path.join(dst,'{}'.format(os.path.basename(pth))+'.npy'), array.astype('float32'))
        if save_tif: tifffile.imsave(os.path.join(dst,'{}'.format(os.path.basename(pth))+'.tif'), array.astype('float32'))
        
        #optional 'save_individual'
        if kwargs['save_individual']:
            im = im*imagescale
            a = np.concatenate((np.max(im, axis=0), np.max(array.astype('uint16'), axis=0)), axis=1)
            b = np.concatenate((np.fliplr(np.rot90(np.max(fix_orientation(im, axes=axes), axis=0),k=3)), np.fliplr(np.rot90(np.max(fix_orientation(array.astype('uint16'), axes=axes), axis=0),k=3))), axis=1)
            plt.figure()
            plt.imshow(np.concatenate((b, a), axis=0), cmap=cmap, alpha=1);  plt.axis('off')
            plt.savefig(os.path.join(dst,'{}'.format(os.path.basename(pth))+'.pdf'), dpi=300, transparent=True)
            plt.close()
        
        #find all nonzero pixels in resampled for elastix volume
        print('   finding nonzero pixels for voxel counts...\n')      
        nz = np.nonzero(array)
        nonzeros.append(zip(*nz)) #<-for pooled image 

        transformfile = transformfiles[i]
        #apply resizing point transform
        txtflnm = point_transform_due_to_resizing(array, chtype = 'injch', **dct)    
        #run transformix on points
        points_file = point_transformix(txtflnm, transformfile)
        #map transformed points to atlas
        tdf = transformed_pnts_to_atlas(points_file, ann, ch_type = 'injch', point_or_index = None, id_table_pth = id_table, **dct) #map to whichever atlas you registered to atlas
        if i == 0: 
            df = tdf.copy()
            countcol = 'count' if 'count' in df.columns else 'cell_count'
            df.drop([countcol], axis=1, inplace=True)
        df[os.path.basename(pth)] = tdf[countcol
         
    #cell counts to csv                           
    df.to_csv(os.path.join(dst,'voxel_counts.csv'))
    print('\n\nCSV file of cell counts, saved as {}\n\n\n'.format(os.path.join(dst,'voxel_counts.csv')))                
        
    return df
    
#%%
def optimize_inj_detect(src, threshold=10, filter_kernel = (5,5,5), dst=False):
    '''Function to test detection parameters
    
    src: path to resized resampled for elastix injection channel volume
    'dst': (optional) path+extension to save image
    
    '''
    if type(src) == str: src = tifffile.imread(src)
    arr = find_site(src, thresh=threshold, filter_kernel=filter_kernel)*45000
    fig = plt.figure()
    fig.add_subplot(1,2,1)
    plt.imshow(np.max(arr, axis=0));  plt.axis('off')
    fig.add_subplot(1,2,2)
    plt.imshow(np.max(src, axis=0), cmap='jet');  plt.axis('off')
    
    if dst: plt.savefig(dst, dpi=300)
    
    return 


def transformed_pnts_to_atlas(points_file, ann, ch_type = 'injch', point_or_index=None, id_table_pth=False, **kwargs):
    '''function to take elastix point transform file and return anatomical locations of those points
    point_or_index=None/point/index: determines which transformix output to use: point is more accurate, index is pixel value(?)
    Elastix uses the xyz convention rather than the zyx numpy convention
    NOTE: this modification does not output out a single excel file, but a data frame
    
    Inputs
    -----------
    points_file = 
    ch_type = 'injch' or 'cellch'
    allen_id_table_pth (optional) pth to allen_id_table
    ann = annotation file
    
    Returns
    -----------
    df = data frame containing voxel counts
    
    '''   
    kwargs = load_kwargs(**kwargs)
    #####inputs 
    assert type(points_file)==str
    
    if point_or_index==None:
        point_or_index = 'OutputPoint'
    elif point_or_index == 'point':
        point_or_index = 'OutputPoint'
    elif point_or_index == 'index':
        point_or_index = 'OutputIndexFixed'

    #
    vols=kwargs['volumes']
    reg_vol=[xx for xx in vols if xx.ch_type == 'regch'][0]

    ####load files
    if not id_table_pth:
        id_table = pd.read_excel('/jukebox/wang/zahra/lightsheet_copy/supp_files/ls_id_table_w_voxelcounts.xlsx') ##use for determining neuroanatomical locations according to allen
    else:
        id_table = pd.read_excel(id_table_pth)
    #ann = ann ###zyx
    with open(points_file, "rb") as f:                
        lines=f.readlines()
        f.close()

    #####populate post-transformed array of contour centers
    sys.stdout.write('\n{} points detected\n\n'.format(len(lines)))
    arr=np.empty((len(lines), 3))    
    for i in range(len(lines)):        
        arr[i,...]=lines[i].split()[lines[i].split().index(point_or_index)+3:lines[i].split().index(point_or_index)+6] #x,y,z
        
    #optional save out of points
    np.save(kwargs['outputdirectory']+'/injection/zyx_voxels.npy', np.asarray([(z,y,x) for x,y,z in arr]))
        
    pnts = transformed_pnts_to_allen_helper_func(arr, ann); pnt_lst=[xx for xx in pnts if xx != 0]
    
    #check to see if any points where found
    if len(pnt_lst)==0:
        raise ValueError('pnt_lst is empty')
    else:
        sys.stdout.write('\nlen of pnt_lst({})\n\n'.format(len(pnt_lst)))
    
    #generate dataframe with column
    df = count_structure_lister(id_table, *pnt_lst) 
    
    return df
    
    
    
