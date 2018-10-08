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
import re


inputs = [
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an1_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an2_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an3_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an4_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an5_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an6_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an7_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an8_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an9_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an10_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an11_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an12_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an13_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an14_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an15_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an16_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an17_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an18_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an19_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an20_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an21_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an22_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an23_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an24_cb_1d3x',
    '/jukebox/wang/Jess/DREADD_cruslateralization/lightsheet/processed/ymaze_cfos_an25_cb_1d3x'
    ]

pth = '/home/wanglab/Downloads/test.pdf'

# export overview file
if __name__ == '__main__':
    check_registration_injection(pth, inputs, cerebellum_only = True, axis = 1)
    # axis: 0 = saggital; 1 = coronal

#%%
def correct_kwargs(src): 
    '''Temporary adjustment to correct kwargs after setting up folders in step 0 locally.
    
    Input: source path of output directory
    '''
    #import kwargs
    kwargs=load_kwargs(src) 
    
    #change packagedirectory (which defaults to my lightsheet_copy for some reason???)
    kwargs['packagedirectory'] = os.path.join(src,'lightsheet')
    
    #change parameterfolder
#    if src[14:18] == 'Jess': #if these are Jess's cerebellum's
#        kwargs['parameterfolder']=src+'/lightsheet/parameterfolder_cb'
#    elif src[24:33] == 'rat-brody': #if these are rat images - need more stringent registration
#        kwargs['parameterfolder']=src+'/lightsheet/parameterfolder_rat'
#    else:
    kwargs['parameterfolder'] = os.path.join(src,'lightsheet/parameterfolder')
    
    #save kwargs
    save_kwargs(src+'/param_dict.p', **kwargs)
    
    #return kwargs to check
    return kwargs


def check_registration_injection(pth, inputs, cerebellum_only = False, axis = 0):
    '''Function to output registered brain images from processed directories.
    Useful to determine 'bad' registration.
    
    Inputs:
        pth = pdf file path
        inputs = Directories containing processed data
    '''
    pdf_pages = PdfPages(pth) #compiles into multiple pdfs
    
    #iterate through inputs
    for src in inputs:
        
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
        
        vol = [xx for xx in dct['volumes'] if xx.ch_type == 'injch'][0] #get injection volume
    
        #read transformed injch image
        print('Reading registered injection channel image\n     {}'.format(outdr))
        im = tifffile.imread(os.path.dirname(vol.ch_to_reg_to_atlas)+'/result.tif')#.astype('uint8')
         
        print('\nPlotting images...\n')
        figs = plt.figure(figsize=(8.27, 11.69))
        #starting to plot figures
        plt.subplot(131)        
        #plot the result and atlas next to each other
        plt.imshow(atl[200], cmap = 'gray'); plt.axis('off'); plt.title('Atlas', fontsize = 10)        
        plt.subplot(132)
        plt.imshow(reg[200], cmap = 'gray'); plt.axis('off'); plt.title('Registered brain', fontsize = 10)
        
        #plot the max intensity zplane for the injection channel
        plt.subplot(133)
        a = np.max(im, axis = axis) # coronal view = 1; saggital view = 0
        plt.imshow(a, cmap = 'plasma', alpha = 1); plt.axis('off'); plt.title('Injection site', fontsize = 10)
        #fix title
        brainname = re.search('(?<=_)(\w+)(?=_1d3x)', vol.brainname)
        
        if cerebellum_only:
            #add title to page
            plt.text(0.5,0.65, brainname.group(0), transform = figs.transFigure, size = 16) #.group(0)
        else:
            #add title to page
            plt.text(0.1,0.70, brainname.group(0), transform = figs.transFigure, size = 16) 
        
        #done with the page
        pdf_pages.savefig(dpi = 300, bbox_inches = 'tight') 
        
    #write PDF document contain composite of all brains
    pdf_pages.close()
    
    print('Saved as {}'.format(pth))
        

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    