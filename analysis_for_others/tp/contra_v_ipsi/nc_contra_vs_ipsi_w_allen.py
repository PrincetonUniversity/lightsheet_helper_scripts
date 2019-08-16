#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:07:37 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl
from tools.registration.register import transformed_pnts_to_allen_helper_func, count_structure_lister
from tools.registration.register import change_transform_parameter_initial_transform
from tools.registration.transform_list_of_points import create_text_file_for_elastix, modify_transform_files
from tools.registration.transform_list_of_points import point_transformix, unpack_pnts
from tools.utils.io import makedir
from skimage.external import tifffile
from tools.analysis.network_analysis import make_structure_objects
from scipy.ndimage.measurements import center_of_mass
 
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/"

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
        raise Exception("Site not detected, try a lower threshold?")
    elif nlab == 1:
        return labelled.astype(bool)
    elif num_sites_to_keep == 1:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        return labelled == np.argmax(sizes)+1
    else:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        vals = [i+1 for i in np.argsort(sizes)[-num_sites_to_keep:][::-1]]
        return np.in1d(labelled, vals).reshape(labelled.shape)
   
#%%

ann_pth = os.path.join(dst, "atlases/sagittal_allen_ann_25um_iso_60um_edge_80um_ventricular_erosion.tif")

#cut annotation file in middle
ann = tifffile.imread(ann_pth)
plt.imshow(ann[300])
z,y,x = ann.shape
#make sure each halves are same dimension as original ann
ann_left = np.zeros_like(ann)
ann_left[:int(z/2), :, :] = ann[:int(z/2), :, :] #cut in the middle in x
ann_right = np.zeros_like(ann)
ann_right[int(z/2):, :, :] = ann[int(z/2):, :, :]
plt.imshow(ann_left[120])

#collect 
#brains should be in this order as they were saved in this order for inj analysis
brains = ["20180409_jg46_bl6_lob6a_04", "20180608_jg75","20170204_tp_bl6_cri_1750r_03", "20180608_jg72", "20180416_jg56_bl6_lob8_04",
 "20170116_tp_bl6_lob45_ml_11", "20180417_jg60_bl6_cri_04", "20180410_jg52_bl6_lob7_05", "20170116_tp_bl6_lob7_1000r_10",
 "20180409_jg44_bl6_lob6a_02", "20180410_jg49_bl6_lob45_02", "20180410_jg48_bl6_lob6a_01", "20180612_jg80", "20180608_jg71",
 "20170212_tp_bl6_crii_1000r_02", "20170115_tp_bl6_lob6a_rpv_03", "20170212_tp_bl6_crii_2000r_03", "20180417_jg58_bl6_sim_02",
 "20170130_tp_bl6_sim_1750r_03", "20170115_tp_bl6_lob6b_ml_04", "20180410_jg50_bl6_lob6b_03", "20170115_tp_bl6_lob6a_1000r_02",
 "20170116_tp_bl6_lob45_500r_12", "20180612_jg77", "20180612_jg76", "20180416_jg55_bl6_lob8_03", "20170115_tp_bl6_lob6a_500r_01",
 "20170130_tp_bl6_sim_rpv_01", "20170204_tp_bl6_cri_1000r_02", "20170212_tp_bl6_crii_250r_01", "20180417_jg61_bl6_crii_05",
 "20170116_tp_bl6_lob7_ml_08", "20180409_jg47_bl6_lob6a_05"]
    
src = "/jukebox/wang/zahra/h129_contra_vs_ipsi/reg_to_allen"

flds = [os.path.join(src, xx) for xx in brains]

#pool brain names and L/R designation into dict
lr_dist = {}

#get inj vol roundabout way
for fld in flds:
    brain = os.path.basename(fld)
    print(brain)
    inj_pth = os.path.join(os.path.join(fld, "inj_to_reg"), "result.tif")
    inj_vol = tifffile.imread(inj_pth)
    z,y,x = inj_vol.shape

    arr = find_site(inj_vol[:, 412:, :])
    #find center of mass
    z_c,y_c,x_c = center_of_mass(arr)
    #take distance from center to arbitrary "midline" (aka half of z axis)
    dist = z_c-(z/2)
    #save to dict 
    lr_dist[brain] = dist
    
    if dist < 0:
        print("brain {} has a left-sided injection\n".format(brain))
    elif dist > 0:
        print("brain {} has a right-sided injection\n".format(brain))
    else:
        print("brain has an injection close to midline so not considering it rn\n")

#%%

#make structures
#FIXME: for some reason the allen table does not work on this, is it ok to use PMA        
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"

structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)
lr_brains = list(lr_dist.keys())
atl_dst = os.path.join(dst, "pma_to_aba"); makedir(atl_dst)
id_table = pd.read_excel(df_pth)
#%%
#------------------------------------------------------------------------------------------------------------------------------
#NOTE THAT ONLY HAVE TO DO THIS ONCE!!!! DO NOT NEED TO DO AGAIN UNLESS DOUBLE CHECKIHG
#transform points to allen atlas space
nc_lst = brains

#get brains that we actually need to get cell counts from
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts/transformed_points"
post_transformed = [os.path.join(src, os.path.join(xx, "transformed_points/posttransformed_zyx_voxels.npy")) for xx in lr_brains]
transformfiles = ["/jukebox/wang/zahra/aba_to_pma/TransformParameters.0.txt",
                  "/jukebox/wang/zahra/aba_to_pma/TransformParameters.1.txt"]

#collect 
for fl in post_transformed:
    arr = np.load(fl)
    #make into transformix-friendly text file
    brain = os.path.basename(os.path.dirname(os.path.dirname(fl)))
    print(brain)
    transformed_dst = os.path.join(atl_dst, brain); makedir(atl_dst)
    pretransform_text_file = create_text_file_for_elastix(arr, transformed_dst)
        
    #copy over elastix files
    trfm_fl = modify_transform_files(transformfiles, transformed_dst) 
    change_transform_parameter_initial_transform(trfm_fl[0], 'NoInitialTransform')
   
    #run transformix on points
    points_file = point_transformix(pretransform_text_file, trfm_fl[-1], transformed_dst)
    
    #convert registered points into structure counts
    converted_points = unpack_pnts(points_file, transformed_dst) 
    
#%%
#------------------------------------------------------------------------------------------------------------------------------    
#NOTE THAT ONLY HAVE TO DO THIS ONCE!!!! DO NOT NEED TO DO AGAIN UNLESS DOUBLE CHECKIHG

def transformed_cells_to_ann(fld, ann, dst, fl_nm):
    """ consolidating to one function bc then no need to copy/paste """
    dct = {}
    
    for fl in fld:
        converted_points = os.path.join(fl, "posttransformed_zyx_voxels.npy")
        print(converted_points)
        point_lst = transformed_pnts_to_allen_helper_func(np.load(converted_points), ann, order = "ZYX")
        df = count_structure_lister(id_table, *point_lst).fillna(0)
        #for some reason duplicating columns, so use this
        nm_cnt = pd.Series(df.cell_count.values, df.name.values).to_dict()
        fl_name = os.path.basename(fl)
        dct[fl_name]= nm_cnt
        
    #unpack
    index = dct[list(dct.keys())[0]].keys()
    columns = dct.keys()
    data = np.asarray([[dct[col][idx] for idx in index] for col in columns])
    df = pd.DataFrame(data.T, columns=columns, index=index)
    
    #save before adding projeny counts at each level
    df.to_pickle(os.path.join(dst, fl_nm))
    
    return os.path.join(dst, fl_nm)

pma2aba_transformed = [os.path.join(atl_dst, xx) for xx in lr_brains]
#collect counts from right side
right = transformed_cells_to_ann(pma2aba_transformed, ann_right, dst, "nc_right_side_no_prog_at_each_level_allen_atl.p")
#collect counts from left side
left = transformed_cells_to_ann(pma2aba_transformed, ann_left, dst, "nc_left_side_no_prog_at_each_level_allen_atl.p")

#%%

def get_cell_n_density_counts(brains, structure, structures, cells_regions, scale_factor = 0.025):
    """ consolidating to one function bc then no need to copy/paste """
    #get cell counts adn densities
    #get densities for all the structures
    df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx", index_col = None)
    df = df.drop(columns = ["Unnamed: 0"])
    df = df.sort_values(by = ["name"])
    
    #make new dict - for all brains
    cells_pooled_regions = {} #for raw counts
    volume_pooled_regions = {} #for density
    
    for brain in brains:    
        #make new dict - this is for EACH BRAIN
        c_pooled_regions = {}
        d_pooled_regions = {}
        
        for soi in structure:
            try:
                soi = [s for s in structures if s.name==soi][0]
                counts = [] #store counts in this list
                volume = [] #store volume in this list
                for k, v in cells_regions[brain].items():
                    if k == soi.name:
                        counts.append(v)
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])#*(scale_factor**3))
                progeny = [str(xx.name) for xx in soi.progeny]
                #now sum up progeny
                if len(progeny) > 0:
                    for progen in progeny:
                        for k, v in cells_regions[brain].items():
                            if k == progen and progen != "Primary somatosensory area, unassigned, layer 4,5,6":
                                counts.append(v)
                                #add to volume list from LUT
                                volume.append(df.loc[df.name == progen, "voxels_in_structure"].values[0])
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
            except:
                for k, v in cells_regions[brain].items():
                    if k == soi:
                        counts.append(v)                    
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])#*(scale_factor**3))                c_pooled_regions[soi] = np.sum(np.asarray(counts))
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
                        
        #add to big dict
        cells_pooled_regions[brain] = c_pooled_regions
        volume_pooled_regions[brain] = d_pooled_regions
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
    cell_counts_per_brain = np.asarray(cell_counts_per_brain)
    
    volume_per_brain = []
    #initialise dummy var
    i = []
    for k,v in volume_pooled_regions.items():
        dct = volume_pooled_regions[k]
        for j,l in dct.items():
            i.append(l)  
        volume_per_brain.append(i)
        #re-initialise for next
        i = []  
    volume_per_brain = np.asarray(volume_per_brain)*(scale_factor**3)
    #calculate denisty
    density_per_brain = np.asarray([xx/volume_per_brain[i] for i, xx in enumerate(cell_counts_per_brain)])
    
    return cell_counts_per_brain, density_per_brain

#making dictionary of cells by region
cells_regions = pckl.load(open(os.path.join(dst, "nc_right_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

nc_areas = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#RIGHT SIDE
cell_counts_per_brain_right, density_per_brain_right = get_cell_n_density_counts(brains, nc_areas, structures, cells_regions)
#LEFT SIDE
cells_regions = pckl.load(open(os.path.join(dst, "nc_left_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      
cell_counts_per_brain_left, density_per_brain_left = get_cell_n_density_counts(brains, nc_areas, structures, cells_regions)


#%%
#preprocessing into contra/ipsi counts per brain, per structure

nc_left_counts = cell_counts_per_brain_left
nc_right_counts = cell_counts_per_brain_right
nc_density_left = density_per_brain_left
nc_density_right = density_per_brain_right

lrv = list(lr_dist.values())
lr_brains = list(lr_dist.keys())

#dct is just for my sanity, so im not mixing up brains
_ccontra = []; _cipsi = []; _dcontra = []; _dipsi = []
for i in range(len(lr_brains)):
    if lrv[i] > 0: #right
        #counts
        _ccontra.append(nc_left_counts[i])
        _cipsi.append(nc_right_counts[i])
        #density
        _dcontra.append(nc_density_left[i])
        _dipsi.append(nc_density_right[i])
    elif lrv[i] < 0: #left
        #counts
        _ccontra.append(nc_right_counts[i])
        _cipsi.append(nc_left_counts[i])
        #density
        _dcontra.append(nc_density_right[i])
        _dipsi.append(nc_density_left[i])


_ccontra = np.asarray(_ccontra).T; _dcontra = np.asarray(_dcontra).T
_cipsi = np.asarray(_cipsi).T; _dipsi = np.asarray(_dipsi).T
_dratio = np.asarray([_dcontra[i]/_dipsi[i] for i in range(len(_dcontra))])
_cratio = np.asarray([_ccontra[i]/_cipsi[i] for i in range(len(_ccontra))])
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
sort_ccontra = _ccontra.T[np.argsort(_dist, axis = 0)]
sort_cipsi = _cipsi.T[np.argsort(_dist, axis = 0)]
sort_cratio = _cratio.T[np.argsort(_dist, axis = 0)]
sort_dcontra = _dcontra.T[np.argsort(_dist, axis = 0)]
sort_dipsi = _dipsi.T[np.argsort(_dist, axis = 0)]
sort_dratio = _dratio.T[np.argsort(_dist, axis = 0)]

sort_inj = _inj[np.argsort(_dist)]   
sort_brains = np.array(lr_brains)[np.argsort(_dist)]

print(sort_dist.shape)
print(sort_cratio.shape)

#group nc regions in smaller, meta regions
#frontal, medial, posterior?
#['Infralimbic area', 'Prelimbic area', 'Anterior cingulate area',
#       'Frontal pole, cerebral cortex', 'Orbital area', 'Gustatory areas',
#       'Agranular insular area', 'Visceral area', 'Somatosensory areas',
#       'Somatomotor areas', 'Retrosplenial area',
#       'Posterior parietal association areas', 'Visual areas',
#       'Temporal association areas', 'Auditory areas', 'Ectorhinal area',
#       'Perirhinal area']
#group_ind = [[0:7], [8:10], [10:]]

grps = np.array(["Frontal\n(IL,PL,ACC,ORB,FRP,\nGU,AI,VISC)" , "Medial\n(MO,SS)", "Posterior\n(RSP,PTL,TE,PERI,ECT)"])
sort_ccontra_pool = np.asarray([[np.sum(xx[0:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in sort_ccontra])
sort_dcontra_pool = np.asarray([[np.sum(xx[0:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in sort_dcontra])
sort_cipsi_pool = np.asarray([[np.sum(xx[0:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in sort_cipsi])
sort_dipsi_pool = np.asarray([[np.sum(xx[0:7]), np.sum(xx[8:10]), np.sum(xx[10:])] for xx in sort_dipsi])
sort_cratio_pool = np.asarray([sort_ccontra_pool[i]/sort_cipsi_pool[i] for i in range(len(sort_ccontra_pool))])
sort_dratio_pool = np.asarray([sort_dcontra_pool[i]/sort_dipsi_pool[i] for i in range(len(sort_dcontra_pool))])
#%%
## display
fig, axes = plt.subplots(ncols = 1, nrows = 5, figsize = (15,5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,1,1,1,0.3]})

#inj fractions
ax = axes[0]

show = np.fliplr(sort_inj).T

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="x-small")

ax = axes[1]
show = sort_dcontra_pool.T
yaxis = grps

vmin = 0
vmax = 500
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
cb.set_label("Density (Cells/$mm^3$)", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 200:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[2]
show = sort_dipsi_pool.T
yaxis = grps

vmin = 0
vmax = 500
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
cb.set_label("Density (Cells/$mm^3$)", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 200:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)


ax = axes[3]
show = sort_dratio_pool.T
yaxis = grps

vmin = 0.7
vmax = 1.5
cmap = plt.cm.Blues
cmap.set_over("navy")
#colormap
bounds = np.linspace(vmin,vmax,5)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
cb.set_label("Density ratio", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col > 1.5:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra/Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[4]
show = np.asarray([sort_dist])
br = lr_brains 

vmin = -100
vmax = 80
cmap = plt.cm.RdBu_r
cmap.set_over('maroon')
cmap.set_under('midnightblue')
#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=2, aspect=10)
cb.set_label("Left to right", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
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

plt.savefig("/home/wanglab/Desktop/density_ratios.pdf", bbox_inches = "tight")

#%%
## display
fig, axes = plt.subplots(ncols = 1, nrows = 5, figsize = (15,5), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,1,1,1,0.3]})

#inj fractions
ax = axes[0]

show = np.fliplr(sort_inj).T

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%d", 
                  shrink=0.9, aspect=5)
cb.set_label("Cell counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="x-small")

ax = axes[1]
show = sort_ccontra_pool.T
yaxis = grps

vmin = 0
vmax = 2000
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.8, aspect=10)
cb.set_label("Cell count", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 200:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[2]
show = sort_cipsi_pool.T
yaxis = grps

vmin = 0
vmax = 2000
cmap = plt.cm.viridis
cmap.set_over("gold")
#colormap
bounds = np.linspace(vmin,vmax,6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%d", shrink=0.8, aspect=10)
cb.set_label("Cell count", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < 200:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)


ax = axes[3]
show = sort_cratio_pool.T
yaxis = grps

vmin = 0.7
vmax = 1.5
cmap = plt.cm.Blues
cmap.set_over("navy")
#colormap
bounds = np.linspace(vmin,vmax,5)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.8, aspect=10)
cb.set_label("Cell count ratio", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col > 1.5:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# aesthetics
ax.set_xticks(np.arange(len(sort_brains))+.5)
sort_brains = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=6, ha="right")
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra/Ipsi", fontsize="small")
ax.yaxis.set_label_coords(-0.15,0.5)

ax = axes[4]
show = np.asarray([sort_dist])
br = lr_brains 

vmin = -100
vmax = 80
cmap = plt.cm.RdBu_r
cmap.set_over('maroon')
cmap.set_under('midnightblue')
#colormap
bounds = np.linspace(vmin,vmax,4)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=2, aspect=10)
cb.set_label("Left to right", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(False)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
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

plt.savefig("/home/wanglab/Desktop/count_ratios.pdf", bbox_inches = "tight")

#basic statistics for these ratios
from scipy.stats import median_absolute_deviation as mad

df = pd.DataFrame()
#decimal to round by
d = 2
df["median_density_ratio"] = np.round(np.median(sort_dratio_pool, axis = 0), d)
df["mean_density_ratio"] = np.round(np.mean(sort_dratio_pool, axis = 0), d)
df["std_density_ratio"] = np.round(np.std(sort_dratio_pool, axis = 0), d)
df["est_std_density_ratio"] = np.round(mad(sort_dratio_pool, axis = 0)/0.6745, d)
df["median_count_ratio"] = np.round(np.median(sort_cratio_pool, axis = 0), d)
df["mean_count_ratio"] = np.round(np.mean(sort_cratio_pool, axis = 0), d)
df["std_count_ratio"] = np.round(np.std(sort_cratio_pool, axis = 0), d)
df["est_std_count_ratio"] = np.round(mad(sort_cratio_pool, axis = 0)/0.6745, d)

df.index = grps

df.to_csv("/home/wanglab/Desktop/ratio_stats.csv")
