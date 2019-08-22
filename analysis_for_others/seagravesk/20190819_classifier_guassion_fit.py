#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:22:36 2019

@author: wanglab
"""

import pickle, numpy as np, matplotlib.pyplot as plt, os
from scipy.optimize import curve_fit
from scipy import asarray as ar

src = "/jukebox/wang/zahra/kelly_cell_detection_analysis"
ecells = os.path.join(src, "edge_cells.p")

rcells = os.path.join(src, "real_cells.p")

rcells = pickle.load(open(rcells, "rb"), encoding = "latin1")
ecells = pickle.load(open(ecells, "rb"), encoding = "latin1")

#change to see distance in pixels around cell center
cutoff = 3
#look at real cells first
#do for all x,y,z
yprof = []
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        yprof.append(cdct[j]["yprofile"])
#make np array
yprof = np.asarray(yprof)
#normalize
norm_yprof = []
for yy in yprof:    
    norm_yprof.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))    
norm_yprof = np.asarray(norm_yprof)

xprof = []
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        xprof.append(cdct[j]["xprofile"])
#make np array
xprof = np.asarray(xprof)
#normalize
norm_xprof = []
for yy in xprof:    
    norm_xprof.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))
norm_xprof = np.asarray(norm_xprof)

zprof = []
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        zprof.append(cdct[j]["zprofile"])
#make np array
zprof = np.asarray(zprof)
#normalize
norm_zprof = []
for yy in zprof:    
    norm_zprof.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))
norm_zprof = np.asarray(norm_zprof)


#do the same for edge cells

#do for all x,y,z
yprof_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    if type(v) is dict:
        yprof_e.append(cdct["yprofile"])
#make np array
yprof_e = np.asarray(yprof_e)
#normalize
norm_yprof_e = []
for yy in yprof_e:    
    norm_yprof_e.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))    
norm_yprof_e = np.asarray(norm_yprof_e)

xprof_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    if type(v) is dict:
        xprof_e.append(cdct["xprofile"])
#make np array
xprof_e = np.asarray(xprof_e)
#normalize
norm_xprof_e = []
for yy in xprof_e:    
    norm_xprof_e.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))
norm_xprof_e = np.asarray(norm_xprof_e)

zprof_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    if type(v) is dict:
        zprof_e.append(cdct["zprofile"])
#make np array
zprof_e = np.asarray(zprof_e)
#normalize
norm_zprof_e = []
for yy in zprof_e:    
    norm_zprof_e.append(np.asarray([(xx-min(yy))/(max(yy)-min(yy)) for xx in yy]))
norm_zprof_e = np.asarray(norm_zprof_e)


#make into sensible numpy arrays to take mean
norm_z = np.ones((norm_zprof.shape[0], 21))*float("nan")
for i in range(len(norm_zprof)):
    norm_z[i, :len(norm_zprof[i])] = norm_zprof[i]
#mean ignoring nans and last 2 values
norm_z_mean = np.nanmean(norm_z[:,:19], axis = 0)    
norm_y = np.ones((norm_yprof.shape[0], 21))*float("nan")
for i in range(len(norm_yprof)):
    norm_y[i, :len(norm_yprof[i])] = norm_yprof[i]
norm_y_mean = np.nanmean(norm_y, axis = 0)        
norm_x = np.ones((norm_xprof.shape[0], 21))*float("nan")
for i in range(len(norm_xprof)):
    norm_x[i, :len(norm_xprof[i])] = norm_xprof[i]
norm_x_mean = np.nanmean(norm_x, axis = 0)        

#%%

for i in range(len(xprof)):
    try:
        prof = xprof[i][5:15]
        x = ar(range(len(prof)))
        y = prof
        
        # weighted arithmetic mean (corrected - check the section below)
        mean = sum(x * y) / sum(y)
        sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
        
        def Gauss(x, a, x0, sigma):
            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
        
        popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma])
        
        plt.plot(x, y, 'b+:', label='data')
        plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
        plt.xlabel('Distance (pixels)')
        plt.ylabel('Intensity')
    except:
        print(i, xprof[i])
plt.show()

#%%

def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def composite_spectrum(x, # data
                       a, b, # linear baseline
                       a1, x01, sigma1, # 1st line
                       a2, x02, sigma2, # 2nd line
                       a3, x03, sigma3 ): # 3rd line
    return (x*a + b + func(x, a1, x01, sigma1)
                    + func(x, a2, x02, sigma2)
                    + func(x, a3, x03, sigma3))

guess = [650, 800, 700, 1000, 800, 1200, 900, 800, 750, 700]

popt, pcov = curve_fit(composite_spectrum, x, y, p0 = guess)
plt.plot(x, composite_spectrum(x, *popt), 'k', label='Total fit')
plt.plot(x, func(x, *popt[-3:])+x*popt[0]+popt[1], c='r', label='Broad component')
FWHM = round(2*np.sqrt(2*np.log(2))*popt[10],4)
plt.axvspan(popt[9]-FWHM/2, popt[9]+FWHM/2, facecolor='g', alpha=0.3, label='FWHM = %s'%(FWHM))
plt.legend(fontsize=10)
plt.show()
