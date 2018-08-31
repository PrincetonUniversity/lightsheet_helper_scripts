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
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def correct_kwargs(src): 
    '''Temporary adjustment to correct kwargs after setting up folders in step 0 locally.
    
    Input: source path of output directory
    NEED TO ADJUST THIS IN THE MAIN PIPELINE?
    '''from matplotlib.backends.backend_pdf import PdfPages
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
    dct = load_kwargs(src)    
    
    #set output directory
    outdr = dct['outputdirectory']
    
    #set atlas file path
    atl = tifffile.imread(dct['AtlasFile'])
    
    #determine elastix output path
    elastix_out = os.path.join(outdr, 'elastix')
        
    #set registration channel file path
    reg = tifffile.imread(os.path.join(elastix_out, 'result.1.tif'))
    
    vol = [xx for xx in dct['volumes'] if xx.ch_type == 'injch'][0]

    im = tifffile.imread(os.path.dirname(vol.ch_to_reg_to_atlas)+'/result.tif')#.astype('uint8')
    
    with PdfPages('/home/wanglab/Desktop/20180831_christina_cerebellum.pdf') as pdf:
        #plot the result and atlas next to each other
        figs = plt.figure(figsize = (8, 6))
        plt.subplot(131)
        plt.imshow(atl[300], cmap = 'gray'); plt.axis('off'); plt.title('Atlas', fontsize = 10)
        
        plt.subplot(132)
        plt.imshow(reg[300], cmap = 'gray'); plt.axis('off'); plt.title('Registered cerebellum', fontsize = 10)
        
        plt.subplot(133)
        a = np.max(im.astype('uint16'), axis = 0)
        plt.imshow(a, cmap = 'plasma', alpha = 1); plt.axis('off'); plt.title('Injection site', fontsize = 10)
        #FIXME: normalise intensity of figure; find stack with highest % of high intensity pixel
    
        pdf.savefig(dpi = 300, bbox_inches = 'tight')
        plt.close()
    
        
    
    