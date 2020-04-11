#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 19:54:05 2020

@author: wanglab
"""

import pandas as pd, matplotlib.pyplot as plt, os, matplotlib as mpl, json
from skimage.external import tifffile
import SimpleITK as sitk, numpy as np, scipy, matplotlib.patches as mpatches, sys
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.registration.allen_structure_json_to_pandas import annotation_location_to_structure
from tools.analysis.network_analysis import make_structure_objects
sys.path.append("/jukebox/wang/zahra/python/ClearMapCluster")
from ClearMap.Analysis.Voxelization import voxelize
import ClearMap.Analysis.Statistics as stat
import ClearMap.Alignment.Resampling as rsp
import ClearMap.IO.IO as io

#formatting LUTs for annotation overlay?
tab20cmap = [plt.cm.tab20(xx) for xx in range(20)]
tab20cmap_nogray = tab20cmap[:14] + tab20cmap[16:]
tab20cmap_nogray = mpl.colors.ListedColormap(tab20cmap_nogray, name = "tab20cmap_nogray")

def get_progeny(dic,parent_structure,progeny_list):
    """ 
    ---PURPOSE---
    Get a list of all progeny of a structure name.
    This is a recursive function which is why progeny_list is an
    argument and is not returned.
    ---INPUT---
    dic                  A dictionary representing the JSON file 
                         which contains the ontology of interest
    parent_structure     The structure
    progeny_list         The list to which this function will 
                         append the progeny structures. 
    """
    if "msg" in list(dic.keys()): dic = dic["msg"][0]
    
    name = dic.get("name")
    children = dic.get("children")
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get("name")
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get("name")
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 

def make_heatmaps(pth, subdir):
    """ 
    makes clearmap style heatmaps 
    NOTE: do not need to do this if you have already succesfully run your images using ClearMap
    """
    #make heatmaps
    vox_params = {"method": "Spherical", "size": (15, 15, 15), "weights": None}
    
    if subdir in os.listdir(pth):
        points = np.load(os.path.join(pth, subdir+"/posttransformed_zyx_voxels.npy"))
        
        #run clearmap style blurring
        vox = voxelize(np.asarray(points), dataSize = (456, 528, 320), **vox_params)
        dst = os.path.join(pth, subdir+"/cells_heatmap.tif")
        tifffile.imsave(dst, vox.astype("int32"))
    else:
        print("no transformed cells")
    return dst

def consolidate_parents_structures_cfos(id_table, ann, namelist, ontology_file, 
                                        verbose=False):
    """
    Function that generates evenly spaced pixels values based on annotation parents

    Removes 0 from list

    Inputs:
        id_table=path to excel file generated from scripts above
        ann = allen annoation file
        namelist=list of structues names, typically parent structures

    Returns:
        -----------
        nann = new array of bitdepth
        list of value+name combinations
    """

    if type(ann) == str: ann = sitk.GetArrayFromImage(sitk.ReadImage(ann))
    #read in id table
    df_ann = pd.read_excel(id_table)
    
    #remove duplicates and null and root
    namelist = list(set(namelist))
    namelist = [xx for xx in namelist if xx != "null" and xx != "root"]
    namelist.sort()

    #setup
    nann = np.zeros(ann.shape).astype("uint8")
    cmap = [xx for xx in np.linspace(1,255, num=len(namelist))]
    
    #open ontology file
    with open(ontology_file) as json_file:
        ontology_dict = json.load(json_file)
        
    #populate
    for i in range(len(namelist)):
        try:
            nm = namelist[i]
            idnum = df_ann.loc[df_ann.name == nm, "id"].values[0]
            # s = [xx for xx in structures if xx.name==nm][0]
            if verbose: print ("{}, {} of {}, value {}".format(nm, i, 
                            len(namelist)-1, cmap[i]))
            nann[np.where(ann==int(idnum))] = cmap[i]
            #find progeny
            progeny = []; get_progeny(ontology_dict, nm, progeny)
            for ii in progeny:
                if ii[3] != "null": nann[np.where(ann==int(ii[3]))] = cmap[i]
        except Exception as e:
            print(nm, e)
    #sitk.Show(sitk.GetImageFromArray(nann))
    #change nann to have NAN where zeros
    nann = nann.astype("float")
    nann[nann == 0] = "nan"

    return nann, zip(cmap[:], namelist)

def make_2D_overlay_of_heatmaps(src, atl_pth, ann_pth, allen_id_table, 
                                save_dst, ontology_file, positive = True, negative = True,
                                no_structures_to_keep = 20, zstep = 40,
                                colorbar_cutoff = 65):
    
    #read in volumes/dataframes
    vol = tifffile.imread(src)
    atl = tifffile.imread(atl_pth)
    ann = tifffile.imread(ann_pth)
    df_ann = pd.read_excel(allen_id_table)
    
    pvol = vol[:,:,:,1]
    nvol = vol[:,:,:,0]
    atl = np.rot90(np.transpose(atl, [1, 0, 2]), axes = (2,1)) #sagittal to coronal
    ann = np.rot90(np.transpose(ann, [1, 0, 2]), axes = (2,1))
    
    #make sure the atlas and annotation file are the same orientation
    assert atl.shape == ann.shape == nvol.shape
    all_struct_iids = list(np.unique(ann.ravel().astype("int64")))
    
    #collect names of parents
    parent_list = []; all_val_struct_iids = []
    for iid in all_struct_iids:
        try:
            parent = str(df_ann.loc[df_ann["id"] == iid, "parent_name"].values[0])
            if not parent == "nan":
                parent_list.append(parent)
                all_val_struct_iids.append(iid)
            else:
                print("\nid %s is not annotated in the ontology, skipping...\n" % iid)
        except:
            print("\nid %s is not annotated in the ontology, skipping...\n" % iid)
    
    #threshold values
    pvol[pvol!=0.0] = 1.0
    nvol[nvol!=0.0] = 1.0
    
    rngs = range(0, ann.shape[0]+zstep, zstep)
    for iii in range(len(rngs)-1):
        rng = (rngs[iii], rngs[iii+1])
    
        #positive
        if positive:
            print(rng, "positive")
            #get highest
            olist = annotation_location_to_structure(allen_id_table, 
                    zip(*np.nonzero(pvol[rng[0]:rng[1]])), ann[rng[0]:rng[1]])
            srt = sorted(olist, key=lambda x: x[0], reverse=True)
            parent_list = [xx[1] for xx in srt]
            #select only subset
            parent_list=parent_list[:no_structures_to_keep]
            nann, lst = consolidate_parents_structures_cfos(allen_id_table, 
                ann[rng[0]:rng[1]], parent_list, verbose=True, structures=structures)
    
            #make fig
            plt.figure(figsize=(15,18))
            ax = plt.subplot(2,1,1)
            fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap="gray", alpha=1)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            ax.set_title("ABA structures")
            mode = scipy.stats.mode(nann, axis=0, nan_policy="omit") #### THIS IS REALLY IMPORTANT
            most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
            ann_mode = mode[0][0]
            masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
            im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=0.8, vmin=0, vmax=255)
            patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
            plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
            ax.set_anchor("W")
    
            #pvals
            ax = plt.subplot(2,1,2)
            ax.set_title("Number of positively correlated voxels in AP dimension")
            #modify colormap
            my_cmap = plt.cm.viridis(np.arange(plt.cm.RdBu.N))
            #my_cmap = plt.cm.RdYlGn(np.arange(plt.cm.RdBu.N))
            my_cmap[:colorbar_cutoff,:4] = 0.0
            my_cmap = mpl.colors.ListedColormap(my_cmap)
            my_cmap.set_under("w")
            #plot
            fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap="gray", alpha=1)
            #plt.imshow(np.max(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.9) #old way
            plt.imshow(np.sum(pvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
            plt.colorbar()
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            ax.set_anchor("W")
            #plt.tight_layout()
            pdst = os.path.join(save_dst, "positive_overlays_zstep{}".format(zstep))
            if not os.path.exists(pdst): os.mkdir(pdst)
            plt.savefig(os.path.join(pdst, "cfos_z{}-{}.pdf".format(rng[0],rng[1])), dpi=300, transparent=True)
            plt.close()
    
        #negative
        if negative:
            print(rng, "negative")
            #get highest
            olist = annotation_location_to_structure(allen_id_table, zip(*np.nonzero(nvol[rng[0]:rng[1]])), ann[rng[0]:rng[1]])
            srt = sorted(olist, key=lambda x: x[0], reverse=True)
            parent_list = [xx[1] for xx in srt]
            #select only subset
            parent_list=parent_list[0:no_structures_to_keep]
            nann, lst = consolidate_parents_structures_cfos(allen_id_table, ann[rng[0]:rng[1]], parent_list, verbose=True, structures=structures)
    
            #make fig
            plt.figure(figsize=(15,18))
            ax = plt.subplot(2,1,1)
            fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap="gray", alpha=1)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            ax.set_title("ABA structures")
            mode = scipy.stats.mode(nann, axis=0, nan_policy="omit") #### THIS IS REALLY IMPORTANT
            most = list(np.unique(mode[1][0].ravel())); most = sorted(most, reverse=True)
            ann_mode = mode[0][0]
            masked_data = np.ma.masked_where(ann_mode < 0.1, ann_mode)
            im = plt.imshow(masked_data, cmap=tab20cmap_nogray, alpha=0.8, vmin=0, vmax=255)
            patches = [mpatches.Patch(color=im.cmap(im.norm(i[0])), label="{}".format(i[1])) for i in lst]
            plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
            ax.set_anchor("W")
    
            #pvals
            ax = plt.subplot(2,1,2)
            ax.set_title("Number of negatively correlated voxels in AP dimension")
            #modify colormap
            import matplotlib as mpl
            my_cmap = plt.cm.plasma(np.arange(plt.cm.RdBu.N))
            my_cmap[:colorbar_cutoff,:4] = 0.0
            my_cmap = mpl.colors.ListedColormap(my_cmap)
            my_cmap.set_under("w")
            #plot
            fig = plt.imshow(np.max(atl[rng[0]:rng[1]], axis=0), cmap="gray", alpha=1)
            plt.imshow(np.sum(nvol[rng[0]:rng[1]], axis=0), cmap=my_cmap, alpha=0.95, vmin=0, vmax=zstep)
            plt.colorbar()
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            ax.set_anchor("W")
            #plt.tight_layout()
            ndst = os.path.join(save_dst, "negative_overlays_zstep{}".format(zstep))
            if not os.path.exists(ndst): os.mkdir(ndst)
            plt.savefig(os.path.join(ndst, "cfos_z{}-{}.pdf".format(rng[0],rng[1])), dpi=300, transparent=True)
            plt.close()
            
    return parent_list

#%%
if __name__ == "__main__":
   
    ########################MAKE P-VALUE MAPS########################
    #make destination directory
    src = "/jukebox/LightSheetData/falkner-mouse/scooter/clearmap_processed"
    pvaldst = "/jukebox/LightSheetData/falkner-mouse/scooter/pooled_analysis/pvalue_maps"
    if not os.path.exists(pvaldst): os.mkdir(pvaldst)
    
    #set up heatmaps per condition
    mm_du_heatmaps = [os.path.join(src, os.path.join(xx, "cells_heatmap.tif")) 
                      for xx in os.listdir(src) if "mm" in xx]
    fm_du_heatmaps = [os.path.join(src, os.path.join(xx, "cells_heatmap.tif")) 
                      for xx in os.listdir(src) if "fm" in xx]
    mf_du_heatmaps = [os.path.join(src, os.path.join(xx, "cells_heatmap.tif")) 
                      for xx in os.listdir(src) if "mf" in xx]
    
    #read all the heatmaps belonging to each group
    mm = stat.readDataGroup(mm_du_heatmaps)
    fm = stat.readDataGroup(fm_du_heatmaps)
    mf = stat.readDataGroup(mf_du_heatmaps)
    
    #find mean and standard deviation of heatmap in each group
    mm_mean = np.mean(mm, axis = 0)
    mm_std = np.std(mm, axis = 0)
        
    fm_mean = np.mean(fm, axis = 0)
    fm_std = np.std(fm, axis = 0)
    
    mf_mean = np.mean(mf, axis = 0)
    mf_std = np.std(mf, axis = 0)
    
    #write mean and standard dev maps to destination
    tifffile.imsave(os.path.join(pvaldst, "mm_mean.tif"), 
                    np.transpose(mm_mean, [1, 0, 2]).astype("float32"))
    tifffile.imsave(os.path.join(pvaldst, "mm_std.tif"), 
                    np.transpose(mm_std, [1, 0, 2]).astype("float32"))
    
    tifffile.imsave(os.path.join(pvaldst, "fm_mean.tif"), 
                    np.transpose(fm_mean, [1, 0, 2]).astype("float32"))
    tifffile.imsave(os.path.join(pvaldst, "fm_std.tif"), 
                    np.transpose(fm_std, [1, 0, 2]).astype("float32"))
    
    tifffile.imsave(os.path.join(pvaldst, "mf_mean.tif"), 
                    np.transpose(mf_mean, [1, 0, 2]).astype("float32"))
    tifffile.imsave(os.path.join(pvaldst, "mf_std.tif"), 
                    np.transpose(mf_std, [1, 0, 2]).astype("float32"))
    
    #Generate the p-values map
    ##########################
    #first comparison: mm vs. fm
    comparison = ["fm_v_mm"]
    cutoff = 0.01 #set p-value cutoff
    #pcutoff: only display pixels below this level of significance
    pvals, psign = stat.tTestVoxelization(mm.astype("float"), 
                    fm.astype("float"), signed = True, pcutoff = cutoff)
    
    #color the p-values according to their sign 
    #(defined by the sign of the difference of the means between the 2 groups)
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    dst = os.path.join(pvaldst, "pvalues_%s.tif" % comparison)
    io.writeData(dst, rsp.sagittalToCoronalData(pvalsc.astype("float32")));
    
    #second comparison: mm vs. mf
    comparison = ["mf_v_mm"]
    pvals, psign = stat.tTestVoxelization(mm.astype("float"), 
                mf.astype("float"), signed = True, pcutoff = cutoff)
    
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    dst = os.path.join(pvaldst, "pvalues_%s.tif" % comparison)
    io.writeData(dst, rsp.sagittalToCoronalData(pvalsc.astype("float32")))
    
    #third comparison: mf vs. fm
    comparison = ["mf_v_fm"]
    pvals, psign = stat.tTestVoxelization(fm.astype("float"), 
                    mf.astype("float"), signed = True, pcutoff = cutoff)
    
    pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);
    dst = os.path.join(pvaldst, "pvalues_%s.tif" % comparison)
    io.writeData(dst, rsp.sagittalToCoronalData(pvalsc.astype("float32")))

#####################END OF SCRIPT THAT MAKES P-VALUE MAPS####################
#%%
    ########################LOOK AT P-VALUE MAPS IN 2D###############################
    #set destination of p value map you want to analyze
#     atl_src = "/jukebox/LightSheetData/falkner-mouse/allen_atlas/"
#     allen_id_table = os.path.join(atl_src, "allen_id_table_w_voxel_counts.xlsx")
#     ann_pth = os.path.join(atl_src, "annotation_2017_25um_sagittal_forDVscans.nrrd")
#     atl_pth = os.path.join(atl_src, "average_template_25_sagittal_forDVscans.tif")
#     ontology_file = os.path.join(atl_src, "allen.json")

#     comparisons = ["fm_v_mm","mf_v_mm","mf_v_fm"]
#     #make 2d overlays
#     for comp in comparisons:
#         flnm = os.path.join(pvaldst, "pvalues_%s.tif" % comp)
#         make_2D_overlay_of_heatmaps(flnm, atl_pth, ann_pth, allen_id_table, 
#                                     pvaldst, ontology_file, positive = True, negative = True)
    
# #%%        
# #get significant structures by order    
# p_val_maps = [
#         listdirfull(os.path.join(dorsal, "observer_v_control"), "tif")[0],
#         listdirfull(os.path.join(dorsal, "demonstrator_v_observer"), "tif")[0], 
#         listdirfull(os.path.join(dorsal, "demonstrator_v_control"), "tif")[0]
#         ]

# atl = tifffile.imread(atl_pth)
# ann = tifffile.imread(ann_pth)

# atl = np.rot90(np.transpose(atl, [1, 0, 2]), axes = (2,1)) #sagittal to coronal
# ann = np.rot90(np.transpose(ann, [1, 0, 2]), axes = (2,1)) #sagittal to coronal

# df = pd.DataFrame()

# for src in p_val_maps:
#     print(src)
    
#     vol = tifffile.imread(src)
#     #positive correlation
#     pvol = vol[:,:,:,1]
    
#     pix_values = list(np.unique(ann.ravel().astype("float64")))
    
#     #collect names of parents OR ##optionally get sublist
#     #make structures to find parents
#     structures = make_structure_objects(allen_id_table)
#     parent_list = list(set([yy.parent[1] for xx in pix_values for yy in structures if xx == yy.idnum]))
    
#     #find counts of highest labeled areas
#     olist = annotation_location_to_structure(allen_id_table, args=zip(*np.nonzero(pvol)), ann=ann)
#     srt = sorted(olist, key=lambda x: x[0], reverse=True)
#     parent_list = [xx[1] for xx in srt]
    
#     if "dorsal" in src:
#         column_name = os.path.basename(src)[8:-4]+"_dorsal_up_positive_corr"
#     else:
#         column_name = os.path.basename(src)[8:-4]+"_ventral_up_positive_corr"
        
#     df[column_name] = pd.Series(parent_list)
#     #negative correlation            
#     nvol = vol[:,:,:,0]
#     pix_values = list(np.unique(ann.ravel().astype("float64")))

#     structures = make_structure_objects(allen_id_table)
#     parent_list = list(set([yy.parent[1] for xx in pix_values for yy in structures if xx == yy.idnum]))
    
#     olist = annotation_location_to_structure(allen_id_table, args=zip(*np.nonzero(nvol)), ann=ann)
#     srt = sorted(olist, key=lambda x: x[0], reverse=True)
#     parent_list = [xx[1] for xx in srt]
    
#     if "dorsal" in src:
#         column_name = os.path.basename(src)[8:-4]+"_dorsal_up_negative_corr"
#     else:
#         column_name = os.path.basename(src)[8:-4]+"_ventral_up_negative_corr"
    
#     df[column_name] = pd.Series(parent_list)

# #save out
# dfdst = os.path.join(pvaldst, "significant_structures_based_on_pvalue_maps.xlsx")       
# df.to_excel(dfdst, index = None)
    