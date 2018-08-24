#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 08:21:27 2018

@author: wanglab
"""
import os
from tools.utils.io import load_kwargs, save_kwargs
import matplotlib.pyplot as plt; plt.ion()
from skimage.external import tifffile


def correct_kwargs(src): 
    '''Temporary adjustment to correct kwargs after setting up folders in step 0 locally.
    
    Input: source path of output directory
    NEED TO ADJUST THIS IN THE MAIN PIPELINE?
    '''
    #import kwargs
    kwargs=load_kwargs(src) 
    
    #change packagedirectory (which defaults to my lightsheet_copy for some reason???)
    kwargs['packagedirectory']=os.path.join(src,'lightsheet')
    
    #change parameterfolder
#    if src[14:18] == 'Jess': #if these are Jess's cerebellum's
#        kwargs['parameterfolder']=src+'/lightsheet/parameterfolder_cb'
#    elif src[24:33] == 'rat-brody': #if these are rat images - need more stringent registration
#        kwargs['parameterfolder']=src+'/lightsheet/parameterfolder_rat'
#    else:
    kwargs['parameterfolder']=os.path.join(src,'lightsheet/parameterfolder')
    
    #save kwargs
    save_kwargs(src+'/param_dict.p', **kwargs)
    
    #return kwargs to check
    return kwargs

def cerebellum_injection(src):
    '''Function to output registered cerebellum images from processed directories, along with injection channel transformix.
    Useful to determine 'bad' registration and fix individual parameters.
    
    Inputs:
        src = Directory containing processed data (in lab bucket)
    '''
    #load kwargs
    kwargs=load_kwargs(src)
    
    #set output directory
    outdr = kwargs['outputdirectory']
    
    #set atlas file path
    atl = kwargs['AtlasFile']
    
    #determine elastix output path
    elastix_out = os.path.join(outdr, 'elastix')
    
    #determine injection transform path
    transform_out=os.path.join(elastix_out, [zz for zz in os.listdir(elastix_out) if zz.endswith('resized_ch01')][0])
    
    #find transformix result tif
    transform_vol = os.path.join(transform_out, [xx for xx in os.listdir(transform_out) if xx == 'result.tif'][0])
    
    #plot the result and atlas next to each other
    plt.figure() 
    plt.subplot(121) #1 = row; 2 = columns; plot 1
    plt.imshow(tifffile.imread(transform_vol)[175], cmap='gray') #FIXME: normalise intensity of figure; find stack with highest % of high intensity pixels
    
    plt.subplot(122)
    plt.imshow(tifffile.imread(atl)[175], cmap='gray')
    
    plt.savefig(os.path.join(outdr,'combined_registered_volumes_inj.pdf'), transparent=True)
    
    
        
    
    