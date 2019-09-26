#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 13:21:07 2019

@author: wanglab
"""

import matplotlib.colors, os, tifffile
import matplotlib.pyplot as mpl, numpy as np
import ClearMap.Visualization.Plot as plt
import ClearMap.IO as io
import ClearMap.ImageProcessing.BackgroundRemoval as bgr
from ClearMap.ImageProcessing.Filter.DoGFilter import filterDoG
from ClearMap.ImageProcessing.MaximaDetection import findExtendedMaxima
from ClearMap.ImageProcessing.MaximaDetection import findCenterOfMaxima
from ClearMap.ImageProcessing.CellSizeDetection import detectCellShape
from ClearMap.ImageProcessing.CellSizeDetection import findCellSize,\
                                                findCellIntensity

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"])

inputs = "/jukebox/wang/pisano/conv_net/annotations/all_better_res/h129/input_files"
filename = os.path.join(inputs, "20170204_tp_bl6_cri_1000r_02_1hfds_647_0010na_25msec_z7d5um_10povlap_ch00_z200-400_y1000-1350_x2050-2400.tif")

data=tifffile.imread(filename)
plt.plotTiling(data, inverse=True, z=(8,20))

dataBGR=bgr.removeBackground(data.astype('float'), size=(25,25), verbose=True)
tifffile.imsave("/home/wanglab/Desktop/bckgrd.tif", dataBGR.astype("uint16"))
plt.plotTiling(dataBGR, inverse=True,z=(10,16))

dataDoG=filterDoG(dataBGR, size=(3,10,10), verbose=True)
tifffile.imsave("/home/wanglab/Desktop/dog.tif", dataDoG.astype("uint16"))
plt.plotTiling(dataDoG*10, inverse=True,z=(10,16))

dataMax=findExtendedMaxima(dataDoG, hMax=None, verbose=True, threshold=15)
tifffile.imsave("/home/wanglab/Desktop/maxima.tif", (dataDoG/dataDoG.max()).astype("float32"))
plt.plotOverlayLabel(  dataDoG/dataDoG.max(), dataMax.astype('int'), z=(10,16))

mpl.imshow(np.max(dataDoG/dataDoG.max(), axis=2), "gist_yarg")
mpl.imshow(np.max(dataMax, axis=2), cmap, alpha=0.3) 
cells=findCenterOfMaxima(data, dataMax)

dataShape = detectCellShape(dataDoG, cells.astype(int), threshold=500)
plt.plotOverlayLabel(dataDoG/dataDoG.max(), dataShape, z=(10,16))
cellSizes=findCellSize(dataShape, maxLabel=cells.shape[0])
cellIntensities=findCellIntensity(dataBGR, dataShape, maxLabel=cells.shape[0])
mpl.figure()
mpl.plot(cellSizes, cellIntensities,'.')
mpl.xlabel('cell size [voxel]')
mpl.ylabel('cell intensity [au]')

cells_map = np.zeros_like(data)
for cell in cells.astype(int):
    z,y,x = cell
    cells_map[z,y,x] = 6000
tifffile.imsave("/home/wanglab/Desktop/cells.tif", cells_map.astype("uint16"))
