#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 09:17:20 2019

@author: wanglab
"""

%matplotlib inline
import numpy as np, pandas as pd, os, sys, shutil, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl
from tools.registration.register import elastix_command_line_call, jacobian_command_line_call, change_interpolation_order, transformix_command_line_call, transformed_pnts_to_allen_helper_func, count_structure_lister
from tools.utils.io import listdirfull, makedir, load_memmap_arr, load_np, listall, load_kwargs
from skimage.external import tifffile
from tools.analysis.network_analysis import make_structure_objects
from scipy.ndimage.measurements import center_of_mass

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#USING 60um edge erosion and 80 um ventricular erosion for NC, as edge seems to be the break bpoint. No real effect for ventricular so will keep the same
ann_pth = "/jukebox/wang/pisano/Python/atlas/stepwise_erosion/annotation_sagittal_atlas_20um_iso_60um_edge_erosion_80um_ventricular_erosion.tif"#"/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif"

#cut annotation file in middle
ann = tifffile.imread(ann_pth)
plt.imshow(ann[300])
#make horizontal
ann_h = np.transpose(ann, [2, 1, 0])
plt.imshow(ann_h[120])
z,y,x = ann_h.shape
ann_h_left = ann_h[:, :, :int(x/2)] #cut in the middle in x
ann_h_right = ann_h[:, :, int(x/2):]
plt.imshow(ann_h_left[120])

ann_sag_left = np.transpose(ann_h_left, [2, 1, 0])
ann_sag_right = np.transpose(ann_h_right, [2, 1, 0])

def add_progeny_counts_at_each_level(df, df_pth, ann_pth):
    """
    """
    #make structures
    from tools.analysis.network_analysis import make_structure_objects
    structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
    
    #make copy of df so not to count things multiple times often
    ddf = df.copy()
    ddf[:] = 0
    
    #now do prog
    for s in structures:
        if len(s.progeny)!=0:
            #break
            prog = [xx.name for xx in s.progeny]
            s_vals = df[df.index==s.name].values
            prog_df_vals = df[df.index.isin(prog)].sum(0).values
            sums = s_vals + prog_df_vals
            ddf[ddf.index==s.name] = sums
        else:
            ddf[ddf.index==s.name] = df[df.index==s.name].values
            
    return ddf        

#now get injection site and automatically designate L/R
def find_site(im, thresh=3, filter_kernel=(3,3,3), num_sites_to_keep=1):
    """Find a connected area of high intensity, using a basic filter + threshold + connected components approach
    
    by: bdeverett

    Parameters
    ----------
    img : np.ndarray
        3D stack in which to find site (technically need not be 3D, so long as filter parameter is adjusted accordingly)
    thresh: float
        threshold for site-of-interest intensity, in number of standard deviations above the mean
    filter_kernel: tuple
        kernel for filtering of image before thresholding
    num_sites_to_keep: int, number of injection sites to keep, useful if multiple distinct sites
    
    Returns
    --------
    bool array of volume where coordinates where detected
    """
    from scipy.ndimage.filters import gaussian_filter as gfilt
    from scipy.ndimage import label
    if type(im) == str: im = tifffile.imread(im)

    filtered = gfilt(im, filter_kernel)
    thresholded = filtered > filtered.mean() + thresh*filtered.std() 
    labelled,nlab = label(thresholded)

    if nlab == 0:
        raise Exception('Site not detected, try a lower threshold?')
    elif nlab == 1:
        return labelled.astype(bool)
    elif num_sites_to_keep == 1:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        return labelled == np.argmax(sizes)+1
    else:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        vals = [i+1 for i in np.argsort(sizes)[-num_sites_to_keep:][::-1]]
        return np.in1d(labelled, vals).reshape(labelled.shape)
    
#make structures
df_pth = "/jukebox/wang/pisano/Python/lightsheet/supp_files/ls_id_table_w_voxelcounts.xlsx"

structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
   
#collect 
id_table = pd.read_excel(df_pth)
#brains should be in this order as they were saved in this order for inj analysis
brains = ['20180409_jg46_bl6_lob6a_04',
 '20180608_jg75',
 '20170204_tp_bl6_cri_1750r_03',
 '20180608_jg72',
 '20180416_jg56_bl6_lob8_04',
 '20170116_tp_bl6_lob45_ml_11',
 '20180417_jg60_bl6_cri_04',
 '20180410_jg52_bl6_lob7_05',
 '20170116_tp_bl6_lob7_1000r_10',
 '20180409_jg44_bl6_lob6a_02',
 '20180410_jg49_bl6_lob45_02',
 '20180410_jg48_bl6_lob6a_01',
 '20180612_jg80',
 '20180608_jg71',
 '20170212_tp_bl6_crii_1000r_02',
 '20170115_tp_bl6_lob6a_rpv_03',
 '20170212_tp_bl6_crii_2000r_03',
 '20180417_jg58_bl6_sim_02',
 '20170130_tp_bl6_sim_1750r_03',
 '20170115_tp_bl6_lob6b_ml_04',
 '20180410_jg50_bl6_lob6b_03',
 '20170115_tp_bl6_lob6a_1000r_02',
 '20170116_tp_bl6_lob45_500r_12',
 '20180612_jg77',
 '20180612_jg76',
 '20180416_jg55_bl6_lob8_03',
 '20170115_tp_bl6_lob6a_500r_01',
 '20170130_tp_bl6_sim_rpv_01',
 '20170204_tp_bl6_cri_1000r_02',
 '20170212_tp_bl6_crii_250r_01',
 '20180417_jg61_bl6_crii_05',
 '20170116_tp_bl6_lob7_ml_08',
 '20180409_jg47_bl6_lob6a_05']
    
src = "/jukebox/wang/pisano/tracing_output/antero_4x"

flds = [os.path.join(src, xx) for xx in brains]

#pool brain names and L/R designation into dict
lr_designation = {}
lr_dist = {}

#get inj vol roundabout way
for fld in flds:
    brain = os.path.basename(fld)
    print(brain)
    fl = os.listdir(os.path.join(fld, "elastix"))
    fl.sort()
    inj_pth = os.path.join(os.path.join(fld, "elastix"), fl[0]+"/result.tif")
    inj_vol = tifffile.imread(inj_pth)
    inj_vol_h = np.transpose(inj_vol, [2, 1, 0])
    z,y,x = inj_vol_h.shape

    #cutting off at 423, same as tom's analysis
    arr = find_site(inj_vol_h[:, 423:, :])
    arr_left = arr[:, :, :int(x/2)]
    arr_right = arr[:, :, int(x/2):]
    #find center of mass
    z_c,y_c,x_c = center_of_mass(arr)
    #take distance from center to arbitrary "midline" (aka half of x axis)
    dist = (x/2)-x_c
    #save to dict 
    lr_dist[os.path.basename(fld)] = dist
    #i just need the left vs right weight, not the actual segment, so will not save that
    left = np.sum(arr_left.astype(int)>0)
    right = np.sum(arr_right.astype(int)>0)
    
    #now im going to assign an arbitrary factor, saying that if the one side of non-zero voxels are 5 times greater than the other, collect
    #this brain and its L/R designation; if there is not a 5 fold difference, drop it from analysis for now
    factor = 1
    if left > right*factor:
        print("brain {} has a left-sided injection\n".format(brain))
        lr_designation[os.path.basename(fld)] = "left"
    elif right > left*factor:
        print("brain {} has a right-sided injection\n".format(brain))
        lr_designation[os.path.basename(fld)] = "right"
    else:
        print("brain has an injection close to midline so not considering it rn\n")


        
#get brains that we actually need to get cell counts from
lr_brains = list(lr_designation.keys())
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts/transformed_points"

post_transformed = [os.path.join(src, os.path.join(xx, "transformed_points/posttransformed_zyx_voxels.npy")) for xx in lr_brains]
dct = {}
for fl in post_transformed:
    arr = np.load(fl)
    point_lst = transformed_pnts_to_allen_helper_func(arr, ann_sag_left, order = "ZYX")
    print(os.path.basename(os.path.dirname(os.path.dirname(fl))))
    df = count_structure_lister(id_table, *point_lst).fillna(0)
    #for some reason duplicating columns, so use this
    nm_cnt = pd.Series(df.cell_count.values, df.name.values).to_dict()
    fl_name = os.path.basename(os.path.dirname(os.path.dirname(fl)))
    dct[fl_name]= nm_cnt
#unpack
index = dct[list(dct.keys())[0]].keys()
columns = dct.keys()
data = np.asarray([[dct[col][idx] for idx in index] for col in columns])
df = pd.DataFrame(data.T, columns=columns, index=index)

#mis mapping seems to be "blank space detected cells"

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi"
#save before adding projeny counts at each level
df.to_pickle(os.path.join(dst, "nc_left_side_no_prog_at_each_level.p"))

#collect 
dct = {}
for fl in post_transformed:
    arr = np.load(fl)
    point_lst = transformed_pnts_to_allen_helper_func(arr, ann_sag_right, order = "ZYX")
    print(os.path.basename(os.path.dirname(os.path.dirname(fl))))
    df = count_structure_lister(id_table, *point_lst).fillna(0)
    #for some reason duplicating columns, so use this
    nm_cnt = pd.Series(df.cell_count.values, df.name.values).to_dict()
    fl_name = os.path.basename(os.path.dirname(os.path.dirname(fl)))
    dct[fl_name]= nm_cnt
    
#unpack
index = dct[list(dct.keys())[0]].keys()
columns = dct.keys()
data = np.asarray([[dct[col][idx] for idx in index] for col in columns])
df = pd.DataFrame(data.T, columns=columns, index=index)

#save before adding projeny counts at each level
df.to_pickle(os.path.join(dst, "nc_right_side_no_prog_at_each_level.p"))

#%%
#RIGHT SIDE
#setup
#making dictionary of cells by region
cells_regions = pckl.load(open(os.path.join(dst, "nc_right_side_no_prog_at_each_level.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

structure = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]
#get cell counts
#make new dict - for all brains
cells_pooled_regions = {} #for raw counts

for brain in lr_brains:    
    #make new dict - this is for EACH BRAIN
    c_pooled_regions = {}
    
    for soi in structure:
        try:
            soi = [s for s in structures if s.name==soi][0]
            counts = [] #store counts in this list
            for k, v in cells_regions[brain].items():
                if k == soi.name:
                    counts.append(v)
            progeny = [str(xx.name) for xx in soi.progeny]
            #now sum up progeny
            if len(progeny) > 0:
                for progen in progeny:
                    for k, v in cells_regions[brain].items():
                        if k == progen:
                            counts.append(v)
            c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
        except:
            for k, v in cells_regions[brain].items():
                if k == soi:
                    counts.append(v)                    
            #add to volume list from LUT
            c_pooled_regions[soi] = np.sum(np.asarray(counts))
                    
    #add to big dict
    cells_pooled_regions[brain] = c_pooled_regions

#making the proper array per brain where regions are removed
cell_counts_per_brain = []

#initialise dummy var
i = []
for k,v in cells_pooled_regions.items():
    dct = cells_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    cell_counts_per_brain.append(i)
    #re-initialise for next
    i = []  
    
cell_counts_per_brain_right = np.asarray(cell_counts_per_brain)

#LEFT SIDE
cells_regions = pckl.load(open(os.path.join(dst, "nc_left_side_no_prog_at_each_level.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

structure = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]
#get cell counts
#make new dict - for all brains
cells_pooled_regions = {} #for raw counts

for brain in lr_brains:    
    #make new dict - this is for EACH BRAIN
    c_pooled_regions = {}
    
    for soi in structure:
        try:
            soi = [s for s in structures if s.name==soi][0]
            counts = [] #store counts in this list
            for k, v in cells_regions[brain].items():
                if k == soi.name:
                    counts.append(v)
            progeny = [str(xx.name) for xx in soi.progeny]
            #now sum up progeny
            if len(progeny) > 0:
                for progen in progeny:
                    for k, v in cells_regions[brain].items():
                        if k == progen:
                            counts.append(v)
            c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
        except:
            for k, v in cells_regions[brain].items():
                if k == soi:
                    counts.append(v)                    
            #add to volume list from LUT
            c_pooled_regions[soi] = np.sum(np.asarray(counts))
                    
    #add to big dict
    cells_pooled_regions[brain] = c_pooled_regions

#making the proper array per brain where regions are removed
cell_counts_per_brain = []

#initialise dummy var
i = []
for k,v in cells_pooled_regions.items():
    dct = cells_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    cell_counts_per_brain.append(i)
    #re-initialise for next
    i = []  
    
cell_counts_per_brain_left = np.asarray(cell_counts_per_brain)

#%%
#preprocessing into contra/ipsi counts per brain, per structure
#different nc counts
#contra = left inj/right counts, etc.
#sum across nc
#nc_left_counts = np.asarray([np.sum(xx) for xx in cell_counts_per_brain_left]).T
#nc_right_counts = np.asarray([np.sum(xx) for xx in cell_counts_per_brain_right]).T

#PICK AN AREA HERE
#area = 8 #infralimbic, prelimbic, etc.
#nc_area = structure[area]
nc_left_counts = cell_counts_per_brain_left
nc_right_counts = cell_counts_per_brain_right

lrv = list(lr_designation.values())
lr_brains = list(lr_designation.keys())

#dct is just for my sanity, so im not mixing up brains
contra = {}; ipsi = {}; _contra = []; _ipsi = []
for i in range(len(lr_brains)):
    if lrv[i] == "right":
        contra[brains[i]] = nc_left_counts[i]; _contra.append(nc_left_counts[i])
        ipsi[brains[i]] = nc_right_counts[i]; _ipsi.append(nc_right_counts[i])
    elif lrv[i] == "left":
        contra[brains[i]] = nc_right_counts[i]; _contra.append(nc_right_counts[i])
        ipsi[brains[i]] = nc_left_counts[i]; _ipsi.append(nc_left_counts[i])
        
_contra = np.asarray(_contra).T
_ipsi = np.asarray(_ipsi).T
_ratio = np.asarray([_contra[i]/_ipsi[i] for i in range(len(_contra))])
#make into one
_dist = np.asarray(list(lr_dist.values()))

#injection site analysis
pth = "/jukebox/wang/zahra/modeling/h129/neocortex/data.p"
data = pckl.load(open(pth, "rb"), encoding = "latin1")

brains = data["brainnames"]
primary_pool = data["primary_pool"]
ak_pool = data["cb_regions_pool"]
inj = data["expr_all_as_frac_of_inj_pool"]
 
_inj = np.asarray([inj[i] for i in range(len(inj)) if brains[i] in lr_brains])
_primary_pool = np.asarray([primary_pool[i] for i in range(len(primary_pool)) if brains[i] in lr_brains])

#sort by distance
sort_dist = np.sort(_dist)
sort_contra = _contra.T[np.argsort(_dist, axis = 0)]
sort_ipsi = _ipsi.T[np.argsort(_dist, axis = 0)]
sort_ratio = _ratio.T[np.argsort(_dist, axis = 0)]
sort_inj = _inj[np.argsort(_dist)]   

print(sort_dist.shape)
print(sort_ratio.shape)

#%%

#plotting ALL 
#formatting
fig, axes = plt.subplots(ncols = 1, nrows = 3, figsize = (20,4), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [3,2,1]})

#inj fractions
ax = axes[0]

show = np.fliplr(sort_inj).T
br = np.asarray(lr_brains)[np.argsort(_dist)]

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="x-small")

#show raw counts    
ax = axes[1]

show = np.flipud(np.asarray(sort_ratio.T[7:9]))
br = lr_brains 

vmin = 0.5
vmax = 1.1
cmap = plt.cm.viridis 
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,3)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.2, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 1:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="small")
        
ax.set_yticks(np.arange(2)+.5)
ax.set_yticklabels(["Somatosensory", "Somatomotor"], fontsize="x-small")
ax.set_ylabel("Contra/Ipsi ratios")

ax = axes[2]
show = np.asarray([sort_dist])
br = lr_brains 

vmin = -100
vmax = 80
cmap = plt.cm.RdBu_r
cmap.set_over('maroon')
cmap.set_under('midnightblue')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,4)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=1.5, aspect=5)
cb.set_label("Left to right", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < -75 or col > 70:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")        

#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(br))+.5)
lbls = np.asarray(br)
ax.set_xticklabels(br, rotation=30, fontsize=5, ha="right")

ax.set_yticks(np.arange(1)+.5)
ax.set_yticklabels(["M-L distance"], fontsize="x-small")

dst = "/home/wanglab/Desktop"
plt.savefig(os.path.join(dst,"contra_ipsi_ss_sm_v2.pdf"), bbox_inches = "tight")

#%%
#plotting ALL 
#formatting
fig, axes = plt.subplots(ncols = 1, nrows = 3, figsize = (20,10), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5,0.5]})

#inj fractions
ax = axes[0]

show = np.fliplr(sort_inj).T
br = np.asarray(lr_brains)[np.argsort(_dist)]

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="x-small")

#show raw counts    
ax = axes[1]

show = np.flipud(np.asarray(sort_ratio.T))
br = lr_brains 

vmin = 0.5
vmax = 1.1
cmap = plt.cm.viridis 
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,3)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.2, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 0.5:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")
        
ax.set_yticks(np.arange(len(structure))+.5)
ax.set_yticklabels(structure, fontsize="x-small")
ax.set_ylabel("Contra/Ipsi ratios")

ax = axes[2]
show = np.asarray([sort_dist])
br = lr_brains 

vmin = -100
vmax = 80
cmap = plt.cm.RdBu_r
cmap.set_over('maroon')
cmap.set_under('midnightblue')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,4)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=1.5, aspect=5)
cb.set_label("Left to right", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < -75 or col > 70:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="x-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="x-small")        

#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(br))+.5)
lbls = np.asarray(br)
ax.set_xticklabels(br, rotation=30, fontsize=5, ha="right")

ax.set_yticks(np.arange(1)+.5)
ax.set_yticklabels(["M-L distance"], fontsize="x-small")

dst = "/home/wanglab/Desktop"
plt.savefig(os.path.join(dst,"contra_ipsi_v2.pdf"), bbox_inches = "tight")
