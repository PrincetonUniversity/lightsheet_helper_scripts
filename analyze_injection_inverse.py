#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:07:15 2018

@author: wanglab
"""

import os, sys, numpy as np
from skimage.external import tifffile
import matplotlib as mpl
from tools.imageprocessing.orientation import fix_orientation
from tools.registration.transform import count_structure_lister, transformed_pnts_to_allen_helper_func
from tools.analysis.analyze_injection import orientation_crop_check, optimize_inj_detect
from tools.registration.register import make_inverse_transform, point_transform_due_to_resizing, point_transformix
from tools.utils.io import load_kwargs, listdirfull, makedir
from natsort import natsorted
from collections import Counter
import SimpleITK as sitk, pandas as pd
import matplotlib.pyplot as plt; plt.ion()

if __name__ == '__main__':
    
    #check if reorientation is necessary
    src = '/jukebox/wang/Jess/lightsheet_output/lawrence/lawrence_an3_crus_iDisco_488_647_025na_1hfds_z10um_250msec/elastix/lawrence_an3_crus_iDisco_488_647_025na_1hfds_z10um_250msec_resized_ch01/result.tif'
    src = orientation_crop_check(src, axes = ('2','1','0'), crop = False) #'[:,390:,:]'
    
    #optimize detection parameters for inj det
    optimize_inj_detect(src, threshold=4, filter_kernel = (4,4,4))

    #run
    #suggestion: save_individual=True,
    #then inspect individual brains, which you can then remove bad brains from list and rerun function
    inputlist = [
            '/jukebox/wang/pisano/tracing_output/antero_4x/20180418_rbdg01',
            '/jukebox/wang/pisano/tracing_output/antero_4x/20180418_rbdg02',
            '/jukebox/wang/pisano/tracing_output/antero_4x/20180418_rbdg03',
            '/jukebox/wang/pisano/tracing_output/antero_4x/20180418_rbdg04',
            '/jukebox/wang/pisano/tracing_output/antero_4x/20180418_rbdg05',
            '/jukebox/wang/pisano/tracing_output/antero_4x/20180418_rbdg06',
            '/jukebox/wang/pisano/tracing_output/antero_4x/20180418_rbdg07'
            ]
    
    kwargs = {'inputlist': inputlist,
              'channel': '01',
              'channel_type': 'injch',
              'filter_kernel': (5,5,5),
              'threshold': 10,
              'num_sites_to_keep': 1,
              'injectionscale': 45000, 
              'imagescale': 1,
              'reorientation': ('2','0','1'),
              'crop': '[:,390:,:]',
              'dst': '/home/wanglab/Downloads/test',
              'save_individual': True, 
              'colormap': 'plasma', 
              'atlas': '/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif',
#              'annotation': '/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif',
#              'id_table': '/jukebox/LightSheetTransfer/atlas/allen_id_table.xlsx'
            }
    
    df = pool_injections_inversetransform(**kwargs)
    
#%%
def pool_injections_inversetransform(**kwargs):
    '''Function to pool several injection sites. 
    Assumes that the basic registration AND inverse transform using elastix has been run. If not, runs inverse transform.
   
    Inputs
    -----------
    kwargs:
      'inputlist': inputlist, #list of folders generated previously from software
      'channel': '01', 
      'channel_type': 'injch',
      'filter_kernel': (3,3,3), #gaussian blur in pixels (if registered to ABA then 1px likely is 25um)
      'threshold': 3 (int, value to use for thresholding, this value represents the number of stand devs above the mean of the gblurred image)
      'num_sites_to_keep': int, number of injection sites to keep, useful if multiple distinct sites
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
      'atlas': '/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif', #whole brain atlas
      'annotation':'/jukebox/LightSheetTransfer/atlas/allenatlas/annotation_sagittal_atlas_20um_iso.tif',
      'id_table': '/jukebox/LightSheetTransfer/atlas/allen_id_table.xlsx',
      
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
#    id_table = kwargs['id_table'] if 'id_table' in kwargs else '/jukebox/temp_wang/pisano/Python/lightsheet/supp_files/allen_id_table.xlsx'
    save_array = kwargs['save_array'] if 'save_array' in kwargs else False
    save_tif = kwargs['save_tif'] if 'save_tif' in kwargs else False
    num_sites_to_keep = kwargs['num_sites_to_keep'] if 'num_sites_to_keep' in kwargs else 1
    nonzeros = []
#    ann = sitk.GetArrayFromImage(sitk.ReadImage(kwargs['annotation']))
#    allen_id_table = pd.read_excel(id_table)
    
    
    for i in range(len(inputlist)): #to iteratre through brains
        pth = inputlist[i] #path of each processed brain
        dct = load_kwargs(inputlist[i]) #load kwargs of brain
        outdr = dct['outputdirectory'] #set output directory of processed brain
        inj_vol = [xx for xx in dct['volumes'] if xx.ch_type == 'injch'][0] #set injection channel volume
        
        print('  loading:\n     {}'.format(pth))
        #run find site function to segment inj site using non-registered resampled for elastix volume - pulled directly from analyze_injection.py
        array = find_site(inj_vol.resampled_for_elastix_vol, thresh=kwargs['threshold'], filter_kernel=kwargs['filter_kernel'], num_sites_to_keep = num_sites_to_keep)*injscale
        if save_array: np.save(os.path.join(dst,'{}'.format(os.path.basename(pth))+'.npy'), array.astype('float32'))
        if save_tif: tifffile.imsave(os.path.join(dst,'{}'.format(os.path.basename(pth))+'.tif'), array.astype('float32'))
        
        #optional 'save_individual'
        if kwargs['save_individual']:
            im = tifffile.imread(inj_vol.resampled_for_elastix_vol)
            im = im*imagescale
            a = np.concatenate((np.max(im, axis=0), np.max(array.astype('uint16'), axis=0)), axis=1)
            b = np.concatenate((np.fliplr(np.rot90(np.max(fix_orientation(im, axes=axes), axis=0),k=3)), np.fliplr(np.rot90(np.max(fix_orientation(array.astype('uint16'), axes=axes), axis=0),k=3))), axis=1)
            plt.figure()
            plt.imshow(np.concatenate((b, a), axis=0), cmap=cmap, alpha=1);  plt.axis('off')
            plt.savefig(os.path.join(dst,'{}'.format(os.path.basename(pth))+'.pdf'), dpi=300, transparent=True)
            plt.close()
        
        #cell counts to csv
        print('   finding nonzero pixels for voxel counts...')      
        nz = np.nonzero(array)
        nonzeros.append(zip(*nz)) #<-for pooled image

        if os.path.exists(os.path.join(outdr, 'elastix_inverse_transform')): #find elastix inverse transform folder in outdr
            svlc = os.path.join(outdr, 'elastix_inverse_transform') 
            svlc = os.path.join(svlc, '{}_{}'.format(inj_vol.ch_type, inj_vol.brainname))
            if os.path.exists(svlc): #if injection inverse transform ran
                atlas2reg2sig = os.path.join(svlc, inj_vol.resampled_for_elastix_vol[inj_vol.resampled_for_elastix_vol.rfind('/')+1:-4]+'_atlas2reg2sig')
                posttransformix = os.path.join(atlas2reg2sig, 'posttransformix')
                points_file = os.path.join(posttransformix, 'outputpoints.txt') #find post transformed points file path
            else: #next few lines pulled directly from register.py - runs injection site inverse!!!!!!!
                transformfile = make_inverse_transform([xx for xx in dct['volumes'] if xx.ch_type == 'injch'][0], cores = 6, **dct)        
                #detect injection site  ##FIXME need to define image and pass in appropriate thresh/filter-kernels
                #inj = [xx for xx in dct['volumes'] if xx.ch_type == 'injch'][0]
                #array = find_site(inj.ch_to_reg_to_atlas+'/result.1.tif', thresh=10, filter_kernel=(5,5,5))         
                #array = find_site(inj.resampled_for_elastix_vol, thresh=10, filter_kernel=(5,5,5)).astype(int)           
                #apply resizing point transform
                txtflnm = point_transform_due_to_resizing(array, chtype = 'injch', **dct)    
                #run transformix on points
                points_file = point_transformix(txtflnm, transformfile)
                
        if os.path.exists(points_file): #if transformed points exist
            tdf = transformed_pnts_to_allen(points_file, ch_type = 'injch', point_or_index = None, **dct)
            if i == 0: 
                df = tdf.copy()
                countcol = 'count' if 'count' in df.columns else 'cell_count'
                df.drop([countcol], axis=1, inplace=True)
            df[os.path.basename(pth)] = tdf[countcol]
        else:
            print 'Points file not found. Make sure param_dict.p is set correctly, and basic registration is completed.'           
         
                               
    df.to_csv(os.path.join(dst,'voxel_counts.csv'))
    print('\n\nCSV file of cell counts, saved as {}\n\n\n'.format(os.path.join(dst,'voxel_counts.csv')))                
                
    #condense nonzero pixels
    nzs = [str(x) for xx in nonzeros for x in xx] #this list has duplicates if two brains had the same voxel w label
    c = Counter(nzs)
    arr = np.zeros(im.shape)
    print('Collecting nonzero pixels for pooled image...')
    tick = 0
    #generating pooled array where voxel value = total number of brains with that voxel as positive
    for k, v in c.iteritems():
        k = [int(xx) for xx in k.replace('(','').replace(')','').split(',')]
        arr[k[0], k[1], k[2]] = int(v)
        tick+=1
        if tick % 50000 == 0: print('   {}'.format(tick))
        
    #load atlas and generate final figure
    print('Generating final figure...')      
    atlas = tifffile.imread(kwargs['atlas'])
    array = fix_orientation(arr, axes=axes)
    
    my_cmap = eval('plt.cm.{}(np.arange(plt.cm.RdBu.N))'.format(cmap))
    my_cmap[:1,:4] = 0.0  
    my_cmap = mpl.colors.ListedColormap(my_cmap)
    my_cmap.set_under('w')
    plt.figure()
    plt.imshow(np.max(atlas, axis=0), cmap='gray')
    plt.imshow(np.max(arr, axis=0), alpha=0.99, cmap=my_cmap); plt.colorbar(); plt.axis('off')
    dpi = int(kwargs['dpi']) if 'dpi' in kwargs else 300
    plt.savefig(os.path.join(dst,'heatmap.pdf'), dpi=dpi, transparent=True);
    plt.close()
    
    print('Saved as {}'.format(os.path.join(dst,'heatmap.pdf')))  
        
    return df
    
    
    
    
    
    
    
    
#%%
def find_site(im, thresh=10, filter_kernel=(5,5,5), num_sites_to_keep=1):
    """Find a connected area of high intensity, using a basic filter + threshold + connected components approach
    
    by: bdeverett
    Parameters
    ----------
    img : np.ndarray
        3D stack in which to find site (technically need not be 3D, so long as filter parameter is adjusted accordingly)
    thresh: float
        threshold for site-of-interest intensity, in number of standard deviations above the mean
    filter_kernel: tuple
        kernel for filtering of image before thresholding
    num_sites_to_keep: int, number of injection sites to keep, useful if multiple distinct sites
    
    Returns
    --------
    bool array of volume where coordinates where detected
    """
    from scipy.ndimage.filters import gaussian_filter as gfilt
    from scipy.ndimage import label
    if type(im) == str: im = tifffile.imread(im)

    filtered = gfilt(im, filter_kernel)
    thresholded = filtered > filtered.mean() + thresh*filtered.std() 
    labelled,nlab = label(thresholded)

    if nlab == 0:
        raise Exception('Site not detected, try a lower threshold?')
    elif nlab == 1:
        return labelled.astype(bool)
    elif num_sites_to_keep == 1:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        return labelled == np.argmax(sizes)+1
    else:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        vals = [i+1 for i in np.argsort(sizes)[-num_sites_to_keep:][::-1]]
	return np.in1d(labelled, vals).reshape(labelled.shape)   
    

def transformed_pnts_to_allen(points_file, ch_type = 'injch', point_or_index=None, allen_id_table_pth=False, **kwargs):
    '''function to take elastix point transform file and return anatomical locations of those points
    point_or_index=None/point/index: determines which transformix output to use: point is more accurate, index is pixel value(?)
    Elastix uses the xyz convention rather than the zyx numpy convention
    
    Inputs
    -----------
    points_file = 
    ch_type = 'injch' or 'cellch'
    allen_id_table_pth (optional) pth to allen_id_table
    
    Returns
    -----------
    excelfl = path to excel file
    
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
    if not allen_id_table_pth:
        allen_id_table=pd.read_excel(os.path.join(reg_vol.packagedirectory, 'supp_files/allen_id_table.xlsx')) ##use for determining neuroanatomical locations according to allen
    else:
        allen_id_table = pd.read_excel(allen_id_table_pth)
    ann=sitk.GetArrayFromImage(sitk.ReadImage(kwargs['annotationfile'])) ###zyx
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
        
    pnts=transformed_pnts_to_allen_helper_func(arr, ann); pnt_lst=[xx for xx in pnts if xx != 0]
    
    #check to see if any points where found
    if len(pnt_lst)==0:
        raise ValueError('pnt_lst is empty')
    else:
        sys.stdout.write('\nlen of pnt_lst({})\n\n'.format(len(pnt_lst)))
    
    #generate dataframe with column
    df = count_structure_lister(allen_id_table, *pnt_lst) 
    
    return df
    
    
    