#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 08:21:27 2018

@author: wanglab
"""

import numpy as np, os, cPickle as pickle
import matplotlib.pyplot as plt; plt.ion()
from skimage.external import tifffile
from matplotlib.backends.backend_pdf import PdfPages
os.chdir('/jukebox/wang/zahra/lightsheet_copy/')
from tools.utils.io import load_kwargs, save_kwargs
os.chdir('/jukebox/wang/zahra/clearmap_cluster_copy/')
from ClearMap.cluster.preprocessing import pth_update

########################################################################RUNS IN PYTHON 2###############################################################
def load_clearmap_kwargs(outdr=None, **kwargs):
    '''simple function to load kwargs given an 'outdr'
    
    Inputs:
    -------------
    outdr: (optional) path to folder generated by package
    kwargs
    '''
    if outdr: kwargs = {}; kwargs = dict([('outputdirectory',outdr)])
    
    with open(pth_update(os.path.join(kwargs['outputdirectory'], 'param_dict.p')), 'rb') as pckl:
        kwargs.update(pickle.load(pckl))
        pckl.close()

    '''
    if update:
        vols = kwargs['volumes']
        [vol.add_brainname(vol.outdr[vol.outdr.rfind('/')+1:]) for vol in vols]
        kwargs['volumes'] = vols
        
        pckloc=os.path.join(outdr, 'param_dict.p'); pckfl=open(pckloc, 'wb'); pickle.dump(kwargs, pckfl); pckfl.close()
    ''' 
    return pth_update(kwargs)

def correct_kwargs(src): 
    '''Temporary adjustment to correct kwargs after setting up folders in step 0 locally.
    
    Input: source path of output directory
    '''
    #import kwargs
    kwargs=load_kwargs(src) 
    
    #change packagedirectory (which defaults to my lightsheet_copy for some reason???)
    kwargs['packagedirectory'] = os.path.join(src,'lightsheet')
    
    kwargs['parameterfolder'] = os.path.join(src,'lightsheet/parameterfolder')
    
    #save kwargs
    save_kwargs(os.path.join(src, 'param_dict.p'), **kwargs)
    
    #return kwargs to check
    return kwargs

def check_clearmap_cfos_output(pth, inputs, axis = 0):
    '''Function to output cfos registered images from clearmap processed directories.
    Useful to determine 'bad' registration or cell count.
    
    Inputs:
        pth = pdf file path
        inputs = Directories containing processed data
    '''
    pdf_pages = PdfPages(pth) #compiles into multiple pdfs
    
    #iterate through inputs
    for src in inputs:
        
        #load kwargs
        dct = load_clearmap_kwargs(src)    
        
        #set output directory
        outdr = dct['outputdirectory']
        
        #determine elastix cfos to auto output path
        elastix_out = os.path.join(outdr, 'clearmap_cluster_output/elastix_auto_to_atlas')
        #set registered volume path
        auto_2_atl = tifffile.imread(os.path.join(elastix_out, 'result.1.tif'))
    
        #read image
        print('Reading autofluo image\n     {}'.format(outdr))
        auto = tifffile.imread(os.path.join(outdr, 'clearmap_cluster_output/autofluo_for_cfos_resampled.tif'))#.astype('uint8')

        print('Reading signal channel image\n     {}'.format(outdr))
        cfos = tifffile.imread(os.path.join(outdr, 'clearmap_cluster_output/cfos_resampled.tif'))#.astype('uint8')
         
        print('\nPlotting images...\n')
        figs = plt.figure(figsize=(8.27, 11.69))
        #starting to plot figures
        plt.subplot(131)        
        #plot the result and atlas next to each other
        plt.imshow(auto[350], cmap = 'gray'); plt.axis('off'); plt.title('Autofluorescent image', fontsize = 10)        
        plt.subplot(132)
        plt.imshow(auto_2_atl[221], cmap = 'gray'); plt.axis('off'); plt.title('Registered Atlas', fontsize = 10)
        plt.subplot(133)
        a = np.max(cfos, axis = axis) # coronal view = 1; saggital view = 0
        plt.imshow(a, cmap = 'plasma', alpha = 1); plt.axis('off'); plt.title('Viral Spread', fontsize = 10)

        #add title to page
        plt.text(0.1,0.70, os.path.basename(src), transform = figs.transFigure, size = 16) 
        
        #done with the page
        pdf_pages.savefig(dpi = 300, bbox_inches = 'tight') 
        
    #write PDF document contain composite of all brains
    pdf_pages.close()
    
    print('Saved as {}'.format(pth))
    
def check_registration_injection(pth, inputs, cerebellum_only = True, axis = 0):
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
        plt.imshow(atl[300], cmap = 'gray'); plt.axis('off'); plt.title('Atlas', fontsize = 10)        
        plt.subplot(132)
        plt.imshow(reg[300], cmap = 'gray'); plt.axis('off'); plt.title('Registered brain', fontsize = 10)
        
        #plot the max intensity zplane for the injection channel
        plt.subplot(133)
        a = np.max(im, axis = axis) # coronal view = 1; saggital view = 0
        plt.imshow(a, cmap = 'plasma', alpha = 1); plt.axis('off'); plt.title('Injection site', fontsize = 10)
        #fix title

#        brainname = re.search('(?<=_)(\w+)(?=_1d3x)', vol.brainname)
        
        if cerebellum_only:
            #add title to page
            plt.text(0.5,0.65, os.path.basename(src), transform = figs.transFigure, size = 16) #.group(0)
        else:
            #add title to page
            plt.text(0.1,0.70, os.path.basename(src), transform = figs.transFigure, size = 16) 
        
        #done with the page
        pdf_pages.savefig(dpi = 300, bbox_inches = 'tight') 
        
    #write PDF document contain composite of all brains
    pdf_pages.close()
    
    print('Saved as {}'.format(pth))
        
#%%

# export overview file
if __name__ == '__main__':
    
#    src = '/jukebox/wang/Jess/lightsheet_output/201812_development/cerebellum'
#
#    inputs = [os.path.join(src, xx) for xx in os.listdir(src) if xx[-5:] == 'crus1']; inputs.sort()
#    
#    pth = '/jukebox/wang/Jess/lightsheet_output/201812_development/201812_development_cerebellum_crus1.pdf'
#
#    check_registration_injection(pth, inputs, cerebellum_only = True, axis = 1)
    # axis: 0 = saggital; 1 = coronal
    src = '/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/processed'

    inputs = [os.path.join(src, xx) for xx in os.listdir(src)]; inputs.sort()
    
    pth = '/jukebox/LightSheetData/pni_viral_vector_core/201902_promoter_exp_6mo/201902_promoter_exp_6mo_overview.pdf'
    
    check_clearmap_cfos_output(pth, inputs, axis = 0)
    # axis: 0 = saggital; 1 = coronal