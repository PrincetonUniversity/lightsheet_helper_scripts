#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 15:40:09 2019

@author: wanglab
"""

import os
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
# matplotlib inline

#collecting data - originally linux file path
src = "/home/wanglab/mounts/LightSheetTransfer/tp/malaz/raw_data"
#src = r"Z:\tp\malaz\raw_data"

brains = [os.path.join(src, xx) for xx in os.listdir(src) if xx != "fluorescence_cohort"]

#initialize intensity profile array
int_profile = []

for brain in brains:

    print(brain)
    #get raw images from left lightsheet
    imgs = [os.path.join(brain, xx) for xx in os.listdir(brain) if "UltraII_raw_RawDataStack[00 x 00]_C00" in xx
            and "UltraII Filter0000" in xx]; imgs.sort()
    
    print("number of files in imgs: {}".format(len(imgs)))
    
    #set range of the "middle'" of the brain
    zrange = [200, 300]
    
    #parse to get to the middle of the brain
    #have to sort to make sure your slices are oriented the same order that you took the image in
    middle = [xx for xx in imgs if int(xx[-30:-27]) > zrange[0] and int(xx[-30:-27]) <= zrange[1]]; middle.sort() 
    
    #TEST
    #visualize one slice
    slc = sitk.GetArrayFromImage(sitk.ReadImage(middle[0]))
#    plt.figure()
#    plt.imshow(slc)
#    
#    print(slc.shape)
#    line = slc[1200, :]
#    plt.figure()
#    plt.plot(line)
    
    #do this for all the slices
    #init line variable
    #doing this across 5 lines
    number_of_lines = 10
    step = 250
    line = np.zeros((number_of_lines, len(middle), slc.shape[1]))
    
    for i in range(len(middle)):
        slc = sitk.GetArrayFromImage(sitk.ReadImage(middle[i]))
        yval = step
        for j in range(number_of_lines):
            line[j, i, :] = slc[yval, :]
            yval += step
    
    #mean across all the lines
    lines = np.mean(line, axis = 0)
    
    #mean across all those slices
    mean_line = lines.mean(axis = 0)
    
    #set condition
    if "dbe-dcm" in brain:
        condition = "dbe-dcm"
    if "dbe+dcm" or "dbe_dcm" in brain:
        condition = "dbe+dcm"
    if "eci+dcm" or "eci_dcm" in brain:
        condition = "eci+dcm"
    if "eci-dcm" or "eci_no_dcm" in brain:
        condition = "eci-dcm"
#    
#    if "dbe_dcm" in brain:
#        condition = "dbe+dcm"
#    elif "eci_dcm" in brain:
#        condition = "eci+dcm"
#    elif "eci_no_dcm" in brain:
#        condition = "eci-dcm"
        
    #save mean line
    int_profile.append([condition, mean_line])


#%%    
#analysis

#overlay all conditions
plt.figure()    
for i in range(len(int_profile)):
    plt.plot(int_profile[i][1])
plt.xlabel("Distance (pixels)")
plt.ylabel("Mean intensity")
    
#overlay by condition
#get data
dbe_no_dcm = []
dbe_dcm = []    
eci_no_dcm = []
eci_dcm = []
for i in range(len(int_profile)):
    if int_profile[i][0] == "dbe-dcm":
        dbe_no_dcm.append(int_profile[i][1])
    elif int_profile[i][0] == "dbe+dcm":
        dbe_dcm.append(int_profile[i][1])
    elif int_profile[i][0] == "eci+dcm":
        eci_dcm.append(int_profile[i][1])
    elif int_profile[i][0] == "eci-dcm":
        eci_no_dcm.append(int_profile[i][1])
        
#take mean per condition
mean_dbe_no_dcm = np.mean(np.asarray(dbe_no_dcm), axis = 0)
mean_dbe_dcm = np.mean(np.asarray(dbe_dcm), axis = 0)
mean_eci_dcm = np.mean(np.asarray(eci_dcm), axis = 0)
mean_eci_no_dcm = np.mean(np.asarray(eci_no_dcm), axis = 0)

plt.figure()
#overlay on each other
plt.plot(mean_dbe_no_dcm[1000:], label = "DBE-DCM")
plt.plot(mean_dbe_dcm[1000:], label = "DBE+DCM")
plt.plot(mean_eci_dcm[1000:], label = "ECI-DCM")
plt.plot(mean_eci_no_dcm[1000:], label = "ECI+DCM")
plt.legend()
plt.xlabel("Distance (pixels)")
plt.ylabel("Mean intensity")
#ADD CODE TO SAVE MATPLOTLIB FIGURE
save_pth = "/home/wanglab/Desktop/test.pdf"
plt.savefig(save_pth, dpi = 300)

#%%
#calculate slopes per animal per condition
from scipy.stats import linregress
#DBE no DCM
dbe_no_dcm_slopes = []
for animal in dbe_no_dcm:
    linreg = linregress(range(len(animal[1000:])), animal[1000:])
    dbe_no_dcm_slopes.append(linreg.slope)

#DBE with DCM
dbe_dcm_slopes = []
for animal in dbe_dcm:
    linreg = linregress(range(len(animal[1000:])), animal[1000:])
    dbe_dcm_slopes.append(linreg.slope)

#ECI with DCM
eci_dcm_slopes = []
for animal in eci_dcm:
    linreg = linregress(range(len(animal[1000:])), animal[1000:])
    eci_dcm_slopes.append(linreg.slope)    
    
#ECI no DCM
eci_no_dcm_slopes = []
for animal in eci_no_dcm:
    linreg = linregress(range(len(animal[1000:])), animal[1000:])
    eci_no_dcm_slopes.append(linreg.slope)    

#stats
#ANOVA ONE WAY
from scipy.stats import f_oneway

f_oneway(eci_dcm_slopes, dbe_dcm_slopes, dbe_no_dcm_slopes)

#T-TEST
from scipy.stats import ttest_ind

#ttest_ind(eci_no_dcm_slopes, dbe_no_dcm_slopes)
ttest_ind(dbe_dcm_slopes, eci_dcm_slopes)

