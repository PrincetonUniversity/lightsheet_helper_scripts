#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 08:21:27 2018

@author: wanglab
"""
from tools.utils.io import load_kwargs, save_kwargs

def correct_kwargs(src): 
    '''Temporary adjustment to correct kwargs after setting up folders in step 0 locally.
    
    Input: source path of output directory
    FIXME: do not include / at the end of src in this funcion
    NEED TO ADJUST THIS IN THE MAIN PIPELINE?
    '''
    #import kwargs
    kwargs=load_kwargs(src) 
    
    #change packagedirectory (which defaults to my lightsheet_copy for some reason???)
    kwargs['packagedirectory']=src+'/lightsheet' #FIX ME: incorporate os.path.join?
    
    #change parameterfolder
    kwargs['parameterfolder']=src+'/lightsheet/parameterfolder'
    
    #save kwargs
    save_kwargs(src+'/param_dict.p', **kwargs)
    
    #return kwargs to check
    return kwargs
