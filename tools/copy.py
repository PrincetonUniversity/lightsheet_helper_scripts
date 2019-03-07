#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 15:21:31 2018

@author: wanglab
"""

import os
import shutil 

def copy(src, dest):
    
    src_files = os.listdir(src)
    src_files.sort()
    
    for file_name in src_files[9:13]:
        fullsizedata = os.path.join(os.path.join(src, file_name), 'full_sizedatafld')
        
        for dct in os.listdir(fullsizedata):            
            if dct[-4:] != '.txt':
                shutil.copytree(os.path.join(fullsizedata, dct), os.path.join(dest, dct))
            
if __name__ == "__main__":
 
#    src = '/home/wanglab/Documents/old_param_dicts'
#    dest = '/home/wanglab/Desktop/test'
    src = '/jukebox/wang/seagravesk/lightsheet/201710_cfos_left_side_only'
    dest = '/jukebox/LightSheetTransfer/kelly/cfos/'
    copy(src, dest)