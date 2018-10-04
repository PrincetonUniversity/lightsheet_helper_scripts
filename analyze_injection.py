#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:14:48 2017

@author: tpisano
"""
import os, numpy as np
from skimage.external import tifffile
import matplotlib as mpl
from tools.imageprocessing.orientation import fix_orientation
from tools.objectdetection.injdetect import find_site
from tools.registration.transform import count_structure_lister, transformed_pnts_to_allen_helper_func
from tools.utils.io import load_kwargs, listdirfull, makedir
from natsort import natsorted
from collections import Counter
import SimpleITK as sitk, pandas as pd
import matplotlib.pyplot as plt; plt.ion()

if __name__ == '__main__':
    
    #check if reorientation is necessary
    src = '/jukebox/wang/Jess/lightsheet_output/christina_dreadds/christina_an2_iDisco_488_647_026na_1hfds_z10um_250msec/elastix/christina_an2_iDisco_488_647_026na_1hfds_z10um_250msec_resized_ch01/result.tif'
    orientation_crop_check(src, axes = ('2','0','1'))
    
    #optimize detection parameters for inj det
    optimize_inj_detect(src, threshold=6, filter_kernel = (5,5,5)) 
    
    #run
    #suggestion: save_individual=True,
    #then inspect individual brains, which you can then remove bad brains from list and rerun function
    inputlist = [
                '/jukebox/wang/zahra/ymaze_cfos_an15_cb'
           #'/jukebox/wang/Jess/lightsheet_output/christina_dreadds/christina_an1_iDisco_488_647_025na_1hfds_z10um_250msec',
#            '/jukebox/wang/Jess/lightsheet_output/christina_dreadds/christina_an2_iDisco_488_647_026na_1hfds_z10um_250msec',
            #'/jukebox/wang/Jess/lightsheet_output/christina_dreadds/christina_an3_iDisco_488_647_025na_1hfds_z10um_250msec',
            #'/jukebox/wang/Jess/lightsheet_output/christina_dreadds/christina_an5_iDisco_488_647_025na_1hfds_z10um_250msec',
            #'/jukebox/wang/Jess/lightsheet_output/christina_dreadds/christina_an6_iDisco_488_647_025na_1hfds_z10um_250msec',
            #'/jukebox/wang/Jess/lightsheet_output/christina_dreadds/christina_an7_iDisco_488_647_025na_1hfds_z10um_250msec',
#            '/jukebox/wang/Jess/lightsheet_output/christina_dreadds/christina_an8_iDisco_488_647_025na_1hfds_z10um_250msec',
#            '/jukebox/wang/Jess/Christina/lightsheet/processed/crusI/christina_an9_iDisco_488_647_025na_1hfds_z10um_250msec',
            #'/jukebox/wang/Jess/Christina/lightsheet/processed/crusI/christina_an10_iDisco_488_647_025na_1hfds_z10um_250msec'
            ]
            

            
    kwargs = {'inputlist': inputlist,
              'channel': '01',
              'channel_type': 'injch',
              'filter_kernel': (5,5,5), #mli(4,4,4)
              'threshold': 4, #mli(4)
              'injectionscale': 10000, 
              'imagescale': 5,
              'reorientation': ('2','0','1'),
              'crop': False,
              'dst': '/home/wanglab/Desktop/test',
              'save_individual': True, 
              'colormap': 'plasma', 
              'atlas': '/jukebox/LightSheetTransfer/atlas/cb_sagittal_atlas_20um_iso.tif',
              'annotation':'/jukebox/LightSheetTransfer/atlas/cb_annotation_sagittal_atlas_20um_iso.tif',
              'id_table': '/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx'
            }              
              
    
    df = pool_injections_for_analysis(**kwargs)
              
#%%
def orientation_crop_check(src, axes = ('0','1','2'), crop = False):
    '''Function to check orientation and cropping. MaxIPs along 0 axis.
    
      'crop': #use to crop volume, values below assume horizontal imaging and sagittal atlas
                False
                cerebellum: '[:,390:,:]'
                caudal midbrain: '[:,300:415,:]'
                midbrain: '[:,215:415,:]'
                thalamus: '[:,215:345,:]'
                anterior cortex: '[:,:250,:]'
    
    '''
    fig = plt.figure()
    plt.axis('off')
    ax = fig.add_subplot(1,2,1)
    im = tifffile.imread(src)
    plt.imshow(np.max(im, axis=0))
    plt.title('Before reorientation')
    
    ax = fig.add_subplot(1,2,2)
    if crop: im = eval('im{}'.format(crop))
    plt.imshow(np.max(fix_orientation(im, axes=axes), axis=0))
    plt.title('After reorientation')
    return




def pool_injections_for_analysis(**kwargs):
    '''Function to pool several injection sites. Assumes that the basic registration using this software has been run.
    
   
    Inputs
    -----------
    kwargs:
      'inputlist': inputlist, #list of folders generated previously from software
      'channel': '01', 
      'channel_type': 'injch',
      'filter_kernel': (3,3,3), #gaussian blur in pixels (if registered to ABA then 1px likely is 25um)
      'threshold': 3 (int, value to use for thresholding, this value represents the number of stand devs above the mean of the gblurred image)
      'injectionscale': 45000, #use to increase intensity of injection site visualizations generated - DOES NOT AFFECT DATA
      'imagescale': 2, #use to increase intensity of background  site visualizations generated - DOES NOT AFFECT DATA
      'reorientation': ('2','0','1'), #use to change image orientation for visualization only
      'crop': #use to crop volume, values below assume horizontal imaging and sagittal atlas
                False
                cerebellum: '[:,390:,:]'
                caudal midbrain: '[:,300:415,:]'
                midbrain: '[:,215:415,:]'
                thalamus: '[:,215:345,:]'
                anterior cortex: '[:,:250,:]'
      
      'dst': '/home/wanglab/Downloads/test', #save location
      'save_individual': True, #optional to save individual images, useful to inspect brains, which you can then remove bad brains from list and rerun function
      'colormap': 'plasma', 
      'atlas': '/home/wanglab/wang/pisano/Python/allenatlas/average_template_25_sagittal_forDVscans.tif',
      'annotation':'/home/wanglab/wang/pisano/Python/allenatlas/annotation_25_ccf2015_forDVscans.nrrd',
      'id_table': '/home/wanglab/temp_wang/pisano/Python/lightsheet/supp_files/allen_id_table.xlsx',
      
      
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
    id_table = kwargs['id_table'] if 'id_table' in kwargs else '/home/wanglab/temp_wang/pisano/Python/lightsheet/supp_files/allen_id_table.xlsx'
    count_threshold = kwargs['count_threshold'] if 'count_threshold' in kwargs else 10
    nonzeros = []
    ann = sitk.GetArrayFromImage(sitk.ReadImage(kwargs['annotation']))
    if kwargs['crop']: ann = eval('ann{}'.format(kwargs['crop']))   
    allen_id_table=pd.read_excel(id_table)
    
    for i in range(len(inputlist)):
        pth = inputlist[i]
        print('Loading: {}'.format(os.path.basename(pth)))
        dct = load_kwargs(pth)
        vol = [xx for xx in dct['volumes'] if xx.ch_type == kwargs['channel_type'] and xx.channel == kwargs['channel']][0]
        #im = tifffile.imread([xx for xx in listdirfull(vol.ch_to_reg_to_atlas) if 'result.{}.'.format(len(listdirfull(vol.parameterfolder))-1) in xx][0])
        if os.path.exists(vol.ch_to_reg_to_atlas):
            im = tifffile.imread(vol.ch_to_reg_to_atlas)#.astype('uint8')
        else:
            im = tifffile.imread(os.path.dirname(vol.ch_to_reg_to_atlas)+'/result.tif')#.astype('uint8')
        if kwargs['crop']: im = eval('im{}'.format(kwargs['crop']))
        
        #segment
        arr = find_site(im, thresh=kwargs['threshold'], filter_kernel=kwargs['filter_kernel'])*injscale
        
        #optional 'save_individual'
        if kwargs['save_individual'] == True:
            im = im*imagescale
            a=np.concatenate((np.max(im, axis=0), np.max(arr.astype('uint16'), axis=0)), axis=1)
            b=np.concatenate((np.fliplr(np.rot90(np.max(fix_orientation(im, axes=axes), axis=0),k=3)), np.fliplr(np.rot90(np.max(fix_orientation(arr.astype('uint16'), axes=axes), axis=0),k=3))), axis=1)
            plt.figure()
            plt.imshow(np.concatenate((b,a), axis=0), cmap=cmap, alpha=1);  plt.axis('off')
            plt.savefig(os.path.join(dst,'{}'.format(os.path.basename(pth))+'.pdf'), dpi=300, transparent=True)
            plt.close()

        #cell counts to csv
        print('   finding nonzero pixels for voxel counts...')      
        nz = np.nonzero(arr)
        nonzeros.append(zip(*nz)) #<-for pooled image
        pos = transformed_pnts_to_allen_helper_func(np.asarray(zip(*[nz[2], nz[1], nz[0]])), ann)
        tdf = count_structure_lister(allen_id_table, *pos)
        if i == 0: 
            df = tdf.copy()
            df.drop(['cell_count'], axis=1, inplace=True)
        df[os.path.basename(pth)] = tdf['cell_count']
        
    df.to_csv(os.path.join(dst,'voxel_counts.csv'))
    print('\n\nCSV file of cell counts, saved as {}\n\n\n'.format(os.path.join(dst,'voxel_counts.csv')))  
            
    #condense nonzero pixels
    nzs = [str(x) for xx in nonzeros for x in xx] #this list has duplicates if two brains had the same voxel w label
    c = Counter(nzs)
    array = np.zeros(im.shape)
    print('Collecting nonzero pixels for pooled image...')
    tick = 0
    #generating pooled array where voxel value = total number of brains with that voxel as positive
    for k,v in c.iteritems():
        k = [int(xx) for xx in k.replace('(','').replace(')','').split(',')]
        array[k[0], k[1], k[2]] = int(v)
        tick+=1
        if tick % 50000 == 0: print('   {}'.format(tick))
        
    #load atlas and generate final figure
    print('Generating final figure...')      
    atlas = tifffile.imread(kwargs['atlas'])
    arr = fix_orientation(array, axes=axes)
    if kwargs['crop']: atlas = eval('atlas{}'.format(kwargs['crop']))
    atlas = fix_orientation(atlas, axes=axes)
    my_cmap = eval('plt.cm.{}(np.arange(plt.cm.RdBu.N))'.format(cmap))
    my_cmap[:1,:4] = 0.0  
    my_cmap = mpl.colors.ListedColormap(my_cmap)
    my_cmap.set_under('w')
    plt.figure()
    plt.imshow(np.max(atlas, axis=0), cmap='gray')
    plt.imshow(np.max(arr, axis=0), alpha=0.99, cmap=my_cmap); plt.colorbar(); plt.axis('off')
    plt.savefig(os.path.join(dst,'heatmap.pdf'), dpi=300, transparent=True);
    plt.close()
    print('Saved as {}'.format(os.path.join(dst,'heatmap.pdf')))  
        
    return df

def optimize_inj_detect(src, threshold=3, filter_kernel = (3,3,3)):
    '''Function to test detection parameters
    '''
    im = tifffile.imread(src)#.astype('uint16')
    arr = find_site(im, thresh=threshold, filter_kernel=filter_kernel)*45000
    fig = plt.figure()
    fig.add_subplot(1,2,1)
    plt.imshow(np.max(arr, axis=0));  plt.axis('off')
    fig.add_subplot(1,2,2)
    plt.imshow(np.max(im, axis=0), cmap='jet');  plt.axis('off')
    return