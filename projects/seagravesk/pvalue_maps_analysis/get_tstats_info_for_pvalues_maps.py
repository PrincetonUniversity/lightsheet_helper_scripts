#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:51:45 2019

@author: wanglab
"""

import os, numpy as np, sys, pandas as pd
from scipy import stats
sys.path.append("/jukebox/wang/zahra/lightsheet_copy")
from tools.utils.io import listdirfull
sys.path.append("/jukebox/wang/zahra/clearmap_cluster_copy")
import ClearMap.Analysis.Statistics as stat

def ttest_stats(group1, group2, dst, cond):
    """t-Test on differences between the individual voxels in group1 and group2, group is a array of voxelizations"""
    
    #make sure destination dir exists
    if not os.path.exists(dst): os.mkdir(dst)
    
    g1 = stat.readDataGroup(group1);  
    g2 = stat.readDataGroup(group2);  
    
    tvals, pvals = stats.ttest_ind(g1, g2, axis = 0, equal_var = False);

    np.save(os.path.join(dst, "tvals_%s.npy" % cond), tvals)
    np.save(os.path.join(dst, "pvals_%s.npy" % cond), pvals)
    
    print("\nsaved t-statistics and p-values corresponding to registered voxels as numpy arrays \n\n.")
    
    pvals_flat = np.ravel(pvals)
    tvals_flat = np.ravel(tvals)
    pvals_flat = pvals_flat[~np.isnan(pvals_flat)]
    tvals_flat = tvals_flat[~np.isnan(tvals_flat)]
    
    pvals_flat.shape
    df = pd.DataFrame()
    df['tvalues'] = tvals_flat
    df['pvalues'] = pvals_flat
    
    df.to_csv(os.path.join(dst, "tstats_no_nans_%s.csv" % cond))
    
    return
#%%    
if __name__ == "__main__":
    
    #2017-2018 dataset
    ann_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif"
    atl_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
    
    pth1 = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos"
    pth2 = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810"
    
    #get dorsal up brains
    control_dorsal_up = [os.path.join(pth1, xx) for xx in os.listdir(pth1) if "mouse" in xx and "_790_" in xx and xx[0] == "1"] #should be 12
    observers_dorsal_up = [os.path.join(pth1, xx) for xx in os.listdir(pth1) if "observer" in xx and "_790_" in xx and xx[0] == "1"] #should be 13
    demonstrator_dorsal_up = [os.path.join(pth1, xx) for xx in os.listdir(pth1) if "demonstrator" in xx and "_790_" in xx and xx[0] == "1"] #should be 12
    
    #get ventral up brains
    control_ventral_up = [os.path.join(pth2, xx) for xx in os.listdir(pth2) if "mouse" in xx and "_790_" in xx]
    control_ventral_up.append("/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_test_20181009/181009_m37079_mouse1_20171014_790_017na_1hfds_z5um_1000msec_11-46-00") #should be 12
    observers_ventral_up = [os.path.join(pth2, xx) for xx in os.listdir(pth2) if "observer" in xx and "_790_" in xx and xx[0] == "1"]
    demonstrator_ventral_up = [os.path.join(pth2, xx) for xx in os.listdir(pth2) if "demonstrator" in xx and "_790_" in xx and xx[0] == "1"] #should be 13
    exclude =  ['/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810/181023_m37081_observer_20171014_790_017na_1hfds_z5um_1000msec_16-20-11',
                '/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos_201810/181017_f37070_observer_20171007_790_017na_1hfds_z5um_1000msec_16-41-12']
    observers_ventral_up = [xx for xx in observers_ventral_up if xx not in exclude] #should be 13
    
    #combine to get all paths
    pths = control_dorsal_up+observers_dorsal_up+demonstrator_dorsal_up+control_ventral_up+observers_ventral_up+demonstrator_ventral_up

    #analysis for dorsal up brains first    
    #make destination directory
    dst = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/pooled_analysis/2017_2018"
    pvaldst = os.path.join(dst, "pvalue_maps")
    if not os.path.exists(dst): os.mkdir(dst)
    if not os.path.exists(pvaldst): os.mkdir(pvaldst)
    if not os.path.exists(pvaldst+"/dorsal_up"): os.mkdir(pvaldst+"/dorsal_up")
    if not os.path.exists(pvaldst+"/ventral_up"): os.mkdir(pvaldst+"/ventral_up")
    
    ctrl_du_heatmaps = [xx+"/cell_region_assignment_99percentile_no_erosion_20190313/cells_heatmap.tif" for xx in control_dorsal_up] 
    obv_du_heatmaps = [xx+"/cell_region_assignment_99percentile_no_erosion_20190313/cells_heatmap.tif" for xx in observers_dorsal_up] 
    dmn_du_heatmaps = [xx+"/cell_region_assignment_99percentile_no_erosion_20190313/cells_heatmap.tif" for xx in demonstrator_dorsal_up] 
    
    c = stat.readDataGroup(ctrl_du_heatmaps)
    o = stat.readDataGroup(obv_du_heatmaps)
    d = stat.readDataGroup(dmn_du_heatmaps)
        
    #comparisons
    conds = ["observer_v_control", "demonstrator_v_observer", "demonstrator_v_control"]
    comparisons = [[o.astype("float"), c.astype("float")], [d.astype("float"), o.astype("float")], [d.astype("float"), c.astype("float")]]
    final_dst = pvaldst+"/dorsal_up"
    
    for i,comp in enumerate(comparisons):
        print(conds[i])
        ttest_stats(comp[0], comp[1], os.path.join(final_dst, conds[i]+"/multiple_comparisons_output"), conds[i])
        
    #%%
    #analysis for ventral up brains
    ctrl_vu_heatmaps = [xx+"/cell_region_assignment_99percentile_no_erosion_20190313/cells_heatmap.tif" for xx in control_ventral_up] 
    obv_vu_heatmaps = [xx+"/cell_region_assignment_99percentile_no_erosion_20190313/cells_heatmap.tif" for xx in observers_ventral_up] 
    dmn_vu_heatmaps = [xx+"/cell_region_assignment_99percentile_no_erosion_20190313/cells_heatmap.tif" for xx in demonstrator_ventral_up] 
    
    c = stat.readDataGroup(ctrl_vu_heatmaps)
    o = stat.readDataGroup(obv_vu_heatmaps)
    d = stat.readDataGroup(dmn_vu_heatmaps)
    
    #comparisons
    conds = ["observer_v_control", "demonstrator_v_observer", "demonstrator_v_control"]
    comparisons = [[o.astype("float"), c.astype("float")], [d.astype("float"), o.astype("float")], [d.astype("float"), c.astype("float")]]
    final_dst = pvaldst+"/ventral_up"
    
    for i,comp in enumerate(comparisons):
        print(conds[i])
        ttest_stats(comp[0], comp[1], os.path.join(final_dst, conds[i]+"/multiple_comparisons_output"), conds[i])

    #%%
    #2019 dataset
    ann_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans_16bit.tif"
    atl_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
    
    pth = "/jukebox/LightSheetTransfer/kelly/201908_cfos"
    subdir = "cell_region_assignment_99percentile_no_erosion_20190909"
    #combine to get all paths
    pths = listdirfull(pth, "647"); pths.sort()
    
    #analysis for dorsal up brains
    #make destination directory
    dst = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/pooled_analysis/2019"
    pvaldst = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/pooled_analysis/2019/pvalue_maps/"
    
    if not os.path.exists(dst): os.mkdir(dst)
    if not os.path.exists(pvaldst): os.mkdir(pvaldst)
    if not os.path.exists(pvaldst+"/dorsal_up"): os.mkdir(pvaldst+"/dorsal_up")
    
    ctrl_du_heatmaps = [os.path.join(pth, os.path.join(xx, subdir+"/cells_heatmap.tif")) for xx in os.listdir(pth) if "647" in xx and "mouse" in xx]; ctrl_du_heatmaps.sort()
    obv_du_heatmaps = [os.path.join(pth, os.path.join(xx, subdir+"/cells_heatmap.tif")) for xx in os.listdir(pth) if "647" in xx and "observ" in xx]; obv_du_heatmaps.sort()
    dmn_du_heatmaps = [os.path.join(pth, os.path.join(xx, subdir+"/cells_heatmap.tif")) for xx in os.listdir(pth) if "647" in xx and "demons" in xx]; dmn_du_heatmaps.sort()
    
    c = stat.readDataGroup(ctrl_du_heatmaps)
    o = stat.readDataGroup(obv_du_heatmaps)
    d = stat.readDataGroup(dmn_du_heatmaps)
    
    #comparisons
    conds = ["observer_v_control", "demonstrator_v_observer", "demonstrator_v_control"]
    comparisons = [[o.astype("float"), c.astype("float")], [d.astype("float"), o.astype("float")], [d.astype("float"), c.astype("float")]]
    final_dst = pvaldst+"/dorsal_up"
    
    for i,comp in enumerate(comparisons):
        print(conds[i])
        ttest_stats(comp[0], comp[1], os.path.join(final_dst, conds[i]+"/multiple_comparisons_output"), conds[i])
    