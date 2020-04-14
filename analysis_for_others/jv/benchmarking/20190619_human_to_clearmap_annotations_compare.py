#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:28:24 2019

@author: wanglab
"""

from __future__ import division
import os, numpy as np
os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from tools.conv_net.functions.bipartite import pairwise_distance_metrics
from tools.utils.io import listdirfull, load_np, makedir, load_dictionary, save_dictionary
from tools.conv_net.input.read_roi import read_roi_zip

if __name__ == "__main__":
    
    #load points dict
    points_dict = load_dictionary("/home/wanglab/Documents/cfos_inputs/filename_points_dictionary.p")   
        
    print(points_dict.keys())
    
    #setup comparison b/w clearmap and annotator(s)? dafina and tiffancy (highest f1 between them)
    src = "/home/wanglab/mounts/wang/Jess/lightsheet_output/201904_ymaze_cfos/clearmap_accuracy_quantification"
    clearmap_dsets = dict([(os.path.basename(xx)[:-13]+".npy", np.load(os.path.join(src, xx))) for xx in os.listdir(src) if xx[-4:] == ".npy"])
    
    #initialise empty vectors
    tps = []; fps = []; fns = []   
    
    #set voxel cutoff value
    cutoff = 5
    
    for i, vol in enumerate(clearmap_dsets.keys()):

        #set ground truth
        clearmap = clearmap_dsets[vol]
        ground_truth = points_dict[vol]
        
        paired,tp,fp,fn = pairwise_distance_metrics(ground_truth, clearmap, cutoff = cutoff) #returns true positive = tp; false positive = fp; false negative = fn
        
        #f1 per dset
        precision = tp/(tp+fp); recall = tp/(tp+fn) #calculating precision and recall
        f1 = 2*( (precision*recall)/(precision+recall) ) #calculating f1 score
        tps.append(tp); fps.append(fp); fns.append(fn) #append matrix to save all values to calculate f1 score
        
    tp = sum(tps); fp = sum(fps); fn = sum(fns) #sum all the elements in the lists
    precision = tp/(tp+fp); recall = tp/(tp+fn) #calculating precision and recall
    f1 = 2*( (precision*recall)/(precision+recall) ) #calculating f1 score
    
    print ("\n   Finished calculating statistics for set params\n\n\nReport:\n***************************\n\
    Cutoff: {} \n\
    F1 score: {}% \n\
    true positives, false positives, false negatives: {} \n\
    precision: {}% \n\
    recall: {}%\n".format(cutoff, round(f1*100, 2), (tp,fp,fn), round(precision*100, 2), round(recall*100, 2)))