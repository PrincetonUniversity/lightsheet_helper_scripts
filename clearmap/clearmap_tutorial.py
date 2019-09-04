#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 13:21:07 2019

@author: wanglab
"""

import matplotlib.colors
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

filename = "/home/wanglab/mounts/wang/zahra/kelly_cell_detection_analysis/comparison_to_clearmap/annotated_volumes/171209_f37104_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_15-46-13_plane500to519_volume1.tif"

data=io.readData(filename)
plt.plotTiling(data, inverse=True, z=(8,20))

dataBGR=bgr.removeBackground(data.astype('float'), size=(7,7), verbose=True)
plt.plotTiling(dataBGR, inverse=True,z=(10,16))

dataDoG=filterDoG(dataBGR, size=(10,10,4), verbose=True)
plt.plotTiling(dataDoG, inverse=True,z=(10,16))

dataMax=findExtendedMaxima(dataDoG, hMax=None, verbose=True, threshold=15)
plt.plotOverlayLabel(  dataDoG/dataDoG.max(), dataMax.astype('int'), z=(10,16))

mpl.imshow(np.max(dataDoG/dataDoG.max(), axis=2), "gist_yarg")
mpl.imshow(np.max(dataMax, axis=2), cmap, alpha=0.3) 
cells=findCenterOfMaxima(data, dataMax)

dataShape = detectCellShape(dataDoG, cells.astype(int), threshold=15)
plt.plotOverlayLabel(dataDoG/dataDoG.max(), dataShape, z=(10,16))
cellSizes=findCellSize(dataShape, maxLabel=cells.shape[0])
cellIntensities=findCellIntensity(dataBGR, dataShape, maxLabel=cells.shape[0])
mpl.figure()
mpl.plot(cellSizes, cellIntensities,'.')
mpl.xlabel('cell size [voxel]')
mpl.ylabel('cell intensity [au]')