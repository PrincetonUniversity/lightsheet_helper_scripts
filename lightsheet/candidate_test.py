#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 14:49:21 2019

@author: wanglab

"""
import os, numpy as np
import tifffile
from collections import Counter
import pandas as pd
from scipy.ndimage.filters import gaussian_filter as gfilt
from scipy.ndimage import label

#brain: 20170308_tp_bl6_lob8_lat_05

def map_injection_to_atlas(img, dst, ann, id_table):
    
    #segment
    seg = find_site(img, thresh=3, filter_kernel=(3,3,3), num_sites_to_keep=1)*45000
    tifffile.imsave(os.path.join(dst, 'brain_inj.tif'), seg.astype('float32'))

    #cell counts to csv
    nz = np.nonzero(seg)
    nonzeros = []
    nonzeros.append(zip(*nz)) #<-for pooled image

    pnt_lst=[]
    arr = np.asarray(list(zip(*[nz[2], nz[1], nz[0]])))
    for i in range(len(arr)):
        pnt=[int(x) for x in arr[i]]
        pnt_lst.append(ann[pnt[2], pnt[1], pnt[0]]) ###find pixel id; arr=XYZ; ann=ZYX

    #make dictionary of pixel id:#num of the id
    cnt = Counter()
    for i in pnt_lst:
        cnt[i]+=1
    
    #generate df + empty column
    id_table = pd.read_excel(id_table) #df=allen_id_table.assign(count= [0]*len(allen_id_table)) #add count columns to df
    tdf = id_table
    tdf['cell_count'] = 0
    
    #populate cell count in dataframe
    for pix_id, count in cnt.items():
        tdf.loc[tdf.id==pix_id, 'cell_count']=count
        
    df = tdf.copy()
    countcol = 'cell_count'
    df.drop([countcol], axis=1, inplace=True)
    df['brain'] = tdf[countcol]
        
    df.to_csv(os.path.join(dst,'voxel_counts.csv'))


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
    
#%%    
if __name__ == '__main__':
    
    img = tifffile.imread("/home/wanglab/Desktop/candidate_test/brain.tif")
    ann = tifffile.imread("/home/wanglab/Desktop/candidate_test/annotations.tif")
    
    id_table = "/home/wanglab/Desktop/candidate_test/id_table.xlsx"
    dst = "/home/wanglab/Desktop/candidate_test"
    map_injection_to_atlas(img, dst, ann, id_table)