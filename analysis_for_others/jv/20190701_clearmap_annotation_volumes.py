#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:33:15 2019

@author: wanglab
"""

import numpy as np, os, tifffile
from tools.conv_net.utils.io import pairwise_distance_metrics, load_dictionary

src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/clearmap_accuracy_quantification/clearmap_results"

sizes = np.linspace(2, 30, 29).astype(int)

flds = [os.path.join(src, xx) for xx in os.listdir(src)]; flds.sort()

n_cells_per_size = []

for sz in flds:
    arrs = [os.path.join(sz, xx) for xx in os.listdir(sz) if xx[-3:] == "npy"]; arrs.sort()
    n_cells_per_size.append(np.array([len(np.load(xx)) for xx in arrs]))

#%%
import pickle

inputs = "/home/wanglab/mounts/wang/Jess/lightsheet_output/201904_ymaze_cfos/clearmap_accuracy_quantification/raw_inputs"
pth = "/home/wanglab/Documents/cfos_inputs/memmap"
if not os.path.exists(pth): os.mkdir(pth)
impths = [os.path.join(inputs, xx) for xx in os.listdir(inputs) if xx[-3:] == "tif"]; impths.sort()
roipths = [os.path.join(inputs, xx) for xx in os.listdir(inputs) if xx[-3:] == "zip"]; roipths.sort()

for impth, roipth in zip(impths, roipths):
    imstack = tifffile.imread(impth)
    
    #dims
    z,y,x = imstack.shape
        
    #init mem mapped array
    arr = np.lib.format.open_memmap(os.path.join(pth, os.path.basename(impth)[:-4]+".npy"), 
                                    dtype = "uint16", mode = "w+", shape = (2,z,y,x))
    
    #fill array:
    arr[0,...] = imstack
    
    #load rois
    if ".zip" in roipth:
        import zipfile
        with zipfile.ZipFile(roipth) as zf:
            rois = zf.namelist()
        #format ZYX, and remove any rois missaved
        rois = [xx for xx in rois if len(xx) == 18] #remove weird missaved ones?
        rois_formatted = list(zip(*[map(int, xx.replace(".roi","").split("-")[0:3]) for xx in rois if len(xx.split("-"))==3]))
    else:
        from tools.conv_net.input.read_roi import read_roi
        with open(roipth, "rb") as fl:
            rois = read_roi(fl)
        rois_formatted = [tuple(xx) for xx in rois]
    
    if len(rois_formatted)==0:
        rois_formatted = zip(*[map(int, xx.replace(".roi","").split("-")[0:3]) for xx in rois if len(xx.split("-"))==4])
    #populate arr; (NOTE: ImageJ has one-based numerics FOR Z but 0 for YX vs np w zero-based numerics for ZYX)
    arr[1,[xx-1 for xx in rois_formatted[0]], rois_formatted[1], rois_formatted[2]] = 255

remove_bad = True
#check all and make dictionary of points (nx3)
file_points_dct = {}
print("Checking all files have valid input and anns...")
for a in [os.path.join(pth, xx) for xx in os.listdir(pth)]:
    arr = np.load(a)
    sh = np.nonzero(arr[0])[0].shape[0]
    pnts = np.nonzero(arr[1])
    shh = pnts[0].shape[0]
    if sh==0 or shh==0:
        print ("File: {} \n  input images nonzeros=={}, annotation nonzeros=={}".format(a,sh,shh))
        if remove_bad:
            os.remove(a)
            print ("removing")
    else:
        file_points_dct[os.path.basename(a)] = zip(*pnts)
        
#save out points
dst = os.path.join(os.path.dirname(pth), "cfos_points_dictionary.p")

with open(dst, 'wb') as fl:    
    pickle.dump(file_points_dct, fl, protocol=pickle.HIGHEST_PROTOCOL)

print("Saved in {}".format(os.path.join(os.path.dirname(pth), "cfos_points_dictionary.p")))    
#%%
    
#points_dict = [(k,list(v)) for k,v in points_dict.items() if len(list(v)) == 0]
#initialise empty vectors
tps_s = []; fps_s = []; fns_s = []   

#set voxel cutoff value
cutoff = 5

for sz in flds:
    arrs = [os.path.join(sz, xx) for xx in os.listdir(sz) if xx[-3:] == "npy"]; arrs.sort()
    
    #reinit    
    tps = []; fps = []; fns = []   
    
    print(sz)
    for i in range(len(arrs)):
        #load points dict        
        points_dict = load_dictionary("/home/wanglab/Documents/cfos_inputs/cfos_points_dictionary.p")   

        #set ground truth
        ground_truth = list(points_dict[os.path.basename(arrs[i][:-13])+".npy"])
        clearmap = np.load(arrs[i])
    
        paired,tp,fp,fn = pairwise_distance_metrics(ground_truth, clearmap, cutoff = cutoff) #returns true positive = tp; false positive = fp; false negative = fn
    
        tps.append(tp); fps.append(fp); fns.append(fn) #append matrix to save all values to calculate f1 score
    
    tps_s.append(tps); fps_s.append(fps); fns_s.append(fns)
    
tp = sum(tps); fp = sum(fps); fn = sum(fns) #sum all the elements in the lists
precision = tp/(tp+fp); recall = tp/(tp+fn) #calculating precision and recall
f1 = 2*( (precision*recall)/(precision+recall) ) #calculating f1 score

#%%

#left y axis
#n_cells_per_size      
n_cells_per_size = np.asarray(n_cells_per_size)
n_cells_per_size_av = np.mean(n_cells_per_size, axis = 1).astype(int)
    
#right y axis (TRUE POSITIVE RATE) OR RECALL (TP/(TP + FN))
#human_ann_propr 
tps_s = np.asarray(tps_s)
fns_s = np.asarray(fns_s)
fps_s = np.asarray(fps_s)
tps_s_av = np.mean(tps_s, axis = 1).astype(int)
tpr = tps_s_av/n_cells_per_size_av
      
#x axis
#cellszthres
cellszthres = sizes

import matplotlib.pyplot as plt

fig = plt.figure()

ax1 = fig.add_subplot(111)
line1 = ax1.plot(n_cells_per_size_av[2:], 'k')
plt.ylabel("Number of cells detected")

# now, the second axes that shares the x-axis with the ax1
ax2 = fig.add_subplot(111, sharex=ax1, frameon=False)
line2 = ax2.plot(tpr[2:], 'r')
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.ylabel("Proportion of true positive cells detected")
plt.xlabel("Cell size threshold (voxels)") 


##all vols
#for im in range(19):
#    # and the first axes using subplot populated with data 
#    ax1 = fig.add_subplot(111)
#    line1 = ax1.plot(n_cells_per_size[1:, im], 'k')
#    plt.ylabel("Number of cells detected")
#
#    # now, the second axes that shares the x-axis with the ax1
#    ax2 = fig.add_subplot(111, sharex=ax1, frameon=False)
#    line2 = ax2.plot(tps_s[1:, im]/n_cells_per_size[1:, im], 'r')
#    ax2.yaxis.tick_right()
#    ax2.yaxis.set_label_position("right")
#    plt.ylabel("Proportion of true positive cells detected")
#    plt.xlabel("Cell size threshold (voxels)") 

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='k', lw=4),
                Line2D([0], [0], color='r', lw=4)]
    
plt.legend(custom_lines, ['Number of cells detected by ClearMap', 'Proportion of detected cells also annotated by user'],
           bbox_to_anchor=(0.97, -0.15))

plt.savefig(os.path.join(src, "true_positive_rate.pdf"), dpi = 300, bbox_inches = "tight")
