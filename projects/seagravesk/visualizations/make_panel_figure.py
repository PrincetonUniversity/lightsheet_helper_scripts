#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 15:18:36 2019

@author: wanglab
"""

import matplotlib.pyplot as plt, tifffile as tif, numpy as np, os

pths = ['/home/wanglab/Desktop/1wkprimary_1wksecondary_f37073_mouse2_20171010_790_maxip_z500-519_left_ls.tif',
        '/home/wanglab/Desktop/1wkprimary_1wksecondary_f62189_mouse1_20190706_647_maxip_z420-439_left_ls.tif',
         '/home/wanglab/Desktop/2wkprimary_1wksecondary_m57201_mouse1_20190830_555_maxip_z360-379_left_ls.tif',
         '/home/wanglab/Desktop/2wkprimary_1wksecondary_m57208_observ_20190910_647_maxip_z540-559_left_ls.tif',
         '/home/wanglab/Desktop/2wkprimary_2wksecondary_m57208_demons_20190910_555_maxip_z365-384_left_ls.tif',
         '/home/wanglab/Desktop/2wkprimary_2wksecondary_m57201_mouse2_20190830_647_maxip_z380-399_left_ls.tif'
         ]

save_dst = "/home/wanglab/Desktop"

brains = [os.path.basename(pth) for pth in pths]
imgs = np.array([tif.imread(pth) for pth in pths])

y,x = imgs[0].shape
bigarr = np.zeros((y*2, x*3)) #2 columns, 3 rows

#fill the array
bigarr[:y,:x] = imgs[0]
bigarr[:y,x:x*2] = imgs[1]
bigarr[:y,x*2:x*3] = imgs[2]
bigarr[y:y*2,:x] = imgs[3]
bigarr[y:y*2,x:x*2] = imgs[4]
bigarr[y:y*2,x*2:x*3] = imgs[5]

tif.imsave("/home/wanglab/Desktop/panel.tif", bigarr.astype("uint16"))
