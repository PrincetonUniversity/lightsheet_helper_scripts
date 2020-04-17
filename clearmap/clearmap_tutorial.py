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
from ClearMap.ImageProcessing.CellSizeDetection import findCellSize,findCellIntensity

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"])

inputs = "/home/wanglab/LightSheetData/rat-brody/processed/201910_tracing/clearmap/"
filename = os.path.join(inputs, "z265_zpln315-340_x4785_y3793_img.tif")

data=tifffile.imread(filename)

dataBGR=bgr.removeBackground(data.astype('float'), size=(25,25), verbose=True)
tifffile.imsave("/home/wanglab/Desktop/bckgrd.tif", dataBGR.astype("uint16"))

dataDoG=filterDoG(dataBGR, size=(3,10,10), verbose=True)
tifffile.imsave("/home/wanglab/Desktop/dog.tif", dataDoG.astype("uint16"))

dataMax=findExtendedMaxima(dataDoG, hMax=None, verbose=True, threshold=15)
tifffile.imsave("/home/wanglab/Desktop/maxima.tif", dataMax.astype('int').astype("uint16"))

mpl.imshow(np.max(dataDoG/dataDoG.max()*1e3, axis=0), "gist_yarg")
mpl.imshow(np.max(dataMax, axis=0), cmap, alpha=0.4)

cells=findCenterOfMaxima(data, dataMax)

dataShape = detectCellShape(dataDoG, cells.astype(int), threshold=20)
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
