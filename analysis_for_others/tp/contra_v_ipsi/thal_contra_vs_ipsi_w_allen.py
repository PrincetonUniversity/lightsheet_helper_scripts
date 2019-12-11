#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:24:45 2019

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
import matplotlib.colors
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop
 
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/"
fig_dst = "/home/wanglab/Desktop"

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
brains = ["20170410_tp_bl6_lob6a_ml_repro_01", "20160823_tp_bl6_cri_500r_02", "20180417_jg59_bl6_cri_03",
"20170207_db_bl6_crii_1300r_02", "20160622_db_bl6_unk_01", "20161205_tp_bl6_sim_750r_03",
"20180410_jg51_bl6_lob6b_04", "20170419_db_bl6_cri_rpv_53hr", "20170116_tp_bl6_lob6b_lpv_07",
"20170411_db_bl6_crii_mid_53hr", "20160822_tp_bl6_crii_1500r_06", "20160920_tp_bl6_lob7_500r_03",
"20170207_db_bl6_crii_rpv_01", "20161205_tp_bl6_sim_250r_02", "20161207_db_bl6_lob6a_500r_53hr",
"20170130_tp_bl6_sim_rlat_05", "20170115_tp_bl6_lob6b_500r_05", "20170419_db_bl6_cri_mid_53hr",
"20161207_db_bl6_lob6a_850r_53hr", "20160622_db_bl6_crii_52hr_01", "20161207_db_bl6_lob6a_50rml_53d5hr",
"20161205_tp_bl6_lob45_1000r_01", "20160801_db_l7_cri_01_mid_64hr"]    

src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/injection_sites"

imgs = [os.path.join(src, xx+".tif.tif") for xx in brains]

#pool brain names and L/R designation into dict
lr_dist = {}
thal_inj_vol = {}

#get inj vol roundabout way
for img in imgs:
    brain = os.path.basename(img)
    print(brain)
    inj_vol = tifffile.imread(img)
    z,y,x = inj_vol.shape
    
    z_c,y_c,x_c = center_of_mass(inj_vol)
    #take distance from center to arbitrary "midline" (aka half of z axis)
    dist = z_c-(z/2)
    #save to dict 
    lr_dist[brain[:-8]] = dist
    thal_inj_vol[brain[:-8]] = np.sum(inj_vol)
    
    if dist < 0:
        print("brain {} has a left-sided injection\n".format(brain))
    elif dist > 0:
        print("brain {} has a right-sided injection\n".format(brain))
    else:
        print("brain has an injection close to midline so not considering it rn\n")


#make structures
#FIXME: for some reason the allen table does not work on this, is it ok to use PMA...    
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ann_df = pd.read_excel(df_pth)

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
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_thalamus/transformed_points"
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
def transformed_cells_to_allen(fld, ann, dst, fl_nm):
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
right = transformed_cells_to_allen(pma2aba_transformed, ann_right, dst, "thal_right_side_no_prog_at_each_level_allen_atl.p")
#collect counts from left side
left = transformed_cells_to_allen(pma2aba_transformed, ann_left, dst, "thal_left_side_no_prog_at_each_level_allen_atl.p")

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
                                volume.append(df.loc[df.name == progen, "voxels_in_structure"].values[0]/2) #divide by 2 since these are half brains!!!
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
            except Exception as e:
                print(e)
                for k, v in cells_regions[brain].items():
                    if k == soi:
                        counts.append(v)                    
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0]/2)#*(scale_factor**3))                c_pooled_regions[soi] = np.sum(np.asarray(counts))
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
    volume_per_brain = np.asarray(volume_per_brain)
    #calculate denisty
    density_per_brain = np.asarray([xx/(volume_per_brain[i]*(scale_factor**3)) for i, xx in enumerate(cell_counts_per_brain)])
    
    return cell_counts_per_brain, density_per_brain, volume_per_brain

#making dictionary of cells by region
cells_regions = pckl.load(open(os.path.join(dst, "thal_right_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      

#nuclei = ["Pons"]
nuclei = ["Pons", "Thalamus, sensory-motor cortex related", "Thalamus, polymodal association cortex related"]
#nuclei = ["Thalamus", "Ventral posteromedial nucleus of the thalamus", "Ventral posterolateral nucleus of the thalamus",
#          "Ventral anterior-lateral complex of the thalamus", "Ventral medial nucleus of the thalamus", "Anteroventral nucleus of thalamus", 
#          "Reticular nucleus of the thalamus", "Ventral part of the lateral geniculate complex", "Mediodorsal nucleus of thalamus",
#          "Submedial nucleus of the thalamus", "Nucleus of reuniens", "Paraventricular nucleus of the thalamus", 
#          "Central lateral nucleus of the thalamus", "Parafascicular nucleus", "Posterior complex of the thalamus",
#          "Lateral dorsal nucleus of thalamus", "Lateral posterior nucleus of the thalamus", "Lateral habenula"]
#nuclei = ["Ventral anterior-lateral complex of the thalamus", "Ventral medial nucleus of the thalamus", 
#          "Ventral posterolateral nucleus of the thalamus", "Ventral posteromedial nucleus of the thalamus",
#          "Subparafascicular nucleus", "Subparafascicular area",
#          "Peripeduncular nucleus", "Medial geniculate complex", "Dorsal part of the lateral geniculate complex",
#          "Lateral posterior nucleus of the thalamus", "Posterior complex of the thalamus",
#          "Posterior limiting nucleus of the thalamus", "Suprageniculate nucleus", "Anteroventral nucleus of thalamus", 
#          "Anteromedial nucleus", "Anterodorsal nucleus", "Interanteromedial nucleus of the thalamus", 
#          "Interanterodorsal nucleus of the thalamus", "Lateral dorsal nucleus of thalamus", 
#          "Intermediodorsal nucleus of the thalamus", "Mediodorsal nucleus of thalamus", "Submedial nucleus of the thalamus", 
#          "Perireunensis nucleus", "Paraventricular nucleus of the thalamus", "Parataenial nucleus", "Nucleus of reuniens", 
#          "Rhomboid nucleus", "Central medial nucleus of the thalamus", "Paracentral nucleus",
#          "Central lateral nucleus of the thalamus", "Parafascicular nucleus", 
#          "Reticular nucleus of the thalamus", "Ventral part of the lateral geniculate complex", "Epithalamus"]

#RIGHT SIDE
cell_counts_per_brain_right, density_per_brain_right, volume_per_brain_right = get_cell_n_density_counts(brains, nuclei, structures, cells_regions)
#LEFT SIDE
cells_regions = pckl.load(open(os.path.join(dst, "thal_left_side_no_prog_at_each_level_allen_atl.p"), "rb"), encoding = "latin1")
cells_regions = cells_regions.to_dict(orient = "dict")      
cell_counts_per_brain_left, density_per_brain_left, volume_per_brain_left = get_cell_n_density_counts(brains, nuclei, structures, cells_regions)


#%%
#preprocessing into contra/ipsi counts per brain, per structure
scale_factor = 0.025
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
_ratio = np.asarray([_ccontra[i]/_cipsi[i] for i in range(len(_ccontra))])
#make into one
_dist = np.asarray(list(lr_dist.values()))

#injection site analysis
data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data.p"
model_data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/model_data_v2.p"
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")
model_data = pckl.load(open(model_data_pth, "rb"), encoding = "latin1")

brains = data["brainnames"]
primary_pool = data["primary_pool"]
ak_pool = data["cb_regions_pool"]
inj = data["expr_all_as_frac_of_inj_pool"]
primary_lob_n = model_data["primary_lob_n"]

_inj = np.asarray([inj[i] for i in range(len(inj)) if brains[i] in lr_brains])
_primary_pool = np.asarray([primary_pool[i] for i in range(len(primary_pool)) if brains[i] in lr_brains])

#sort by distance
sort_dist = np.sort(_dist)
sort_ccontra = _ccontra.T[np.argsort(_dist, axis = 0)]
sort_cipsi = _cipsi.T[np.argsort(_dist, axis = 0)]
sort_ratio = _ratio.T[np.argsort(_dist, axis = 0)]
sort_dcontra = _dcontra.T[np.argsort(_dist, axis = 0)]
sort_dipsi = _dipsi.T[np.argsort(_dist, axis = 0)]
sort_vox_per_region = volume_per_brain_left[np.argsort(_dist, axis = 0)]
sort_inj = _inj[np.argsort(_dist)]   
sort_brains = np.array(lr_brains)[np.argsort(_dist)]

print(sort_dist.shape)
print(sort_ratio.shape)

#%%
#filter results by contra/ipsi ratio > 1 in thalamus

threshold = 1.5
filter_dcontra = np.array([struct[_ratio[0]>threshold] for struct in _dcontra])
filter_dipsi = np.array([struct[_ratio[0]>threshold] for struct in _dipsi])
filter_ccontra = np.array([struct[_ratio[0]>threshold] for struct in _ccontra])
filter_cipsi = np.array([struct[_ratio[0]>threshold] for struct in _cipsi])

filter_primary_pool = _primary_pool[_ratio[0]>threshold]
filter_ak_pool = ak_pool[np.unique(filter_primary_pool)]
filter_primary_lob_n = np.asarray([np.where(filter_primary_pool == i)[0].shape[0] for i in np.unique(filter_primary_pool)])
#show contra, ipsi, and contra+ipsi heatmaps side by side

#mean percent counts
_pccontra = np.nan_to_num(np.array([xx[1:]/xx[0] for xx in filter_ccontra.T])*100) #note that the whole thalamus is the first element in nthe array
mean_contra_pcounts = np.array([np.mean(_pccontra[np.where(filter_primary_pool == idx)[0]], axis=0) for idx in np.unique(filter_primary_pool)])

#get acronyms of nuclei
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
#choose whether to annotate the numbers in the heatmap
annotate = False
#set label coords
xaxis_label_x,xaxis_label_y = 0.7, 0.1
#set range of colormap
vmin = 0
vmax = 10

fig, axes = plt.subplots(ncols = 3, nrows = 1, figsize = (8,6), sharey = True, gridspec_kw = {"wspace":0, "hspace":0})

ax = axes[0]
show = mean_contra_pcounts.T 

cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,((vmax-vmin)/2)+1)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.1, aspect=10)
cb.set_label("% of thalamic cells", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
# exact value annotations
if annotate:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            if col < 1:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
            else:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        

# xticks
ax.set_xticks(np.arange(len(filter_ak_pool))+.5)
lbls = np.asarray(filter_ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, filter_primary_lob_n)], rotation=30, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="small")
#make x label
ax.set_xlabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(xaxis_label_x,xaxis_label_y)

#ipsi side
ax = axes[1]

_pcipsi = np.array([xx[1:]/xx[0] for xx in filter_cipsi.T])*100
mean_ipsi_pcounts = np.asarray([np.mean(_pcipsi[np.where(filter_primary_pool == idx)[0]], axis=0) for idx in np.unique(filter_primary_pool)])

show = mean_ipsi_pcounts.T 

#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.1, aspect=10)
cb.set_label("% of thalamic cells", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(False)
# exact value annotations
if annotate:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            if col < 1:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
            else:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
# xticks
ax.set_xticks(np.arange(len(filter_ak_pool))+.5)
lbls = np.asarray(filter_ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, filter_primary_lob_n)], rotation=30, fontsize=5, ha="right")
#make x label
ax.set_xlabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(xaxis_label_x,xaxis_label_y)

#combined sides
#ipsi side
ax = axes[2]

_pcounts = np.array([xx[1:]/xx[0] for xx in (filter_cipsi+filter_ccontra).T])*100
mean_pcounts = np.asarray([np.mean(_pcounts[np.where(filter_primary_pool == idx)[0]], axis=0) for idx in np.unique(filter_primary_pool)])

show = mean_pcounts.T 

#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
#cb = pl.colorbar(pc, ax=ax, label="Weight / SE", shrink=0.5, aspect=10)
#cb = pl.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%1i", shrink=0.5, aspect=10)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.1, aspect=10)
cb.set_label("% of thalamic cells", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
if annotate:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            if col < 1:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
            else:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
# xticks
ax.set_xticks(np.arange(len(filter_ak_pool))+.5)
lbls = np.asarray(filter_ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, filter_primary_lob_n)], rotation=30, fontsize=5, ha="right")
#make x label
ax.set_xlabel("Bilateral", fontsize="small")
ax.yaxis.set_label_coords(xaxis_label_x,xaxis_label_y)

plt.savefig(os.path.join(fig_dst, "thal_contra_n_ipsi_mean_pcounts_threshold_%s.pdf" % threshold), bbox_inches = "tight")

plt.close()

#filter results by contra/ipsi ratio > 1 in thalamus

#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
pcounts_descending_order = np.sort(_pccontra)
order = np.argsort(np.mean(_pccontra, axis = 0))
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
sois_descending_pcounts = np.array(short_nuclei)[order]

#boxplots of percent counts - contra
fig, axes = plt.subplots(ncols = 3, nrows = 1, figsize = (15,8), sharey = False, sharex = True, gridspec_kw = {"wspace":0, "hspace":0})

#first boxplot
ax = axes[0]

ax.boxplot(pcounts_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(pcounts_descending_order.T)
for i in range(ngroup):
    ax.scatter(pcounts_descending_order[:,i], 
                y=np.ones(len(pcounts_descending_order[:,0]))*i+1, color = "k", s = 10)


ax.set_xlabel("% of thalamic cells\nContra")
ax.set_ylabel("Thalamic nuclei")

#second boxplot, ipsi
ax = axes[1]

pcounts_descending_order = np.sort(_pcipsi)
order = np.argsort(np.mean(_pcipsi, axis = 0))
sois_descending_pcounts = np.array(short_nuclei)[order]

ax.boxplot(pcounts_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(pcounts_descending_order.T)
for i in range(ngroup):
    ax.scatter(pcounts_descending_order[:,i], 
                y=np.ones(len(pcounts_descending_order[:,0]))*i+1, color = "k", s = 10)


ax.set_xlabel("Ipsi")

#final boxplot, combined
ax = axes[2]

pcounts_descending_order = np.sort(_pcounts)
order = np.argsort(np.mean(_pcounts, axis = 0))
sois_descending_pcounts = np.array(short_nuclei)[order]

ax.boxplot(pcounts_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(pcounts_descending_order.T)
for i in range(ngroup):
    ax.scatter(pcounts_descending_order[:,i], 
                y=np.ones(len(pcounts_descending_order[:,0]))*i+1, color = "k", s = 10)

#label which boxplot belongs to which side
ax.set_xlabel("Bilateral")

plt.savefig(os.path.join(fig_dst, "thal_contra_n_ipsi_pcounts_boxplots_threshold_%s.pdf" % threshold), bbox_inches = "tight")

plt.close()

#do a rank correlation between nuclei for contra vs. ipsi side

import statsmodels.api as sm

ak_vh = np.array(["Vermis", "Hemisphere"])
func = lambda xx: 0 if xx < 2 else 1
filter_primary_pool_vh = np.array([func(xx) for xx in filter_primary_pool])

_pccontra_vermis = _pccontra[np.where(filter_primary_pool_vh == 0)]
_pccontra_hem = _pccontra[np.where(filter_primary_pool_vh == 1)]
_pcipsi_vermis = _pcipsi[np.where(filter_primary_pool_vh == 0)]
_pcipsi_hem = _pcipsi[np.where(filter_primary_pool_vh == 1)]

#vermis
contra_pcounts_descending_order = np.sort(_pccontra_vermis)
contra_order = np.argsort(np.mean(_pccontra_vermis, axis = 0))
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
contra_sois_descending_pcounts = np.array(short_nuclei)[contra_order]
contra_sois_ascending_pcounts = contra_sois_descending_pcounts[::-1]

ipsi_pcounts_descending_order = np.sort(_pcipsi_vermis)
ipsi_order = np.argsort(np.mean(_pcipsi_vermis, axis = 0))
ipsi_sois_descending_pcounts = np.array(short_nuclei)[ipsi_order]
ipsi_sois_ascending_pcounts = ipsi_sois_descending_pcounts[::-1]

contra_ranks = [i for i,nuc in enumerate(contra_sois_ascending_pcounts)]
ipsi_ranks = [i for j,nu in enumerate(ipsi_sois_ascending_pcounts) 
    for i,nuc in enumerate(contra_sois_ascending_pcounts) if ipsi_sois_ascending_pcounts[j] == nuc]

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#size of scatter
size = 70

Y = contra_ranks
X = ipsi_ranks

results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
ax.scatter(y = Y, x = X, s = size)

#plot fit line
ax.plot(mean_slope*range(len(X))+mean_intercept, '--k')    
    
ax.set_xlabel("Contralateral thalamus rank order (vermis injections)")
ax.set_ylabel("Ipsilateral thalamus rank order (vermis injections)")

#make text box
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)


plt.savefig(os.path.join(fig_dst, "thal_contra_v_ipsi_rank_order_vermis_threshold_%s.pdf" % threshold), bbox_inches = "tight")

plt.close()

#hemisphere
contra_pcounts_descending_order = np.sort(_pccontra_hem)
contra_order = np.argsort(np.mean(_pccontra_hem, axis = 0))
short_nuclei = [ann_df.loc[ann_df.name == nuc, "acronym"].values[0] for nuc in nuclei][1:]
contra_sois_descending_pcounts = np.array(short_nuclei)[contra_order]
contra_sois_ascending_pcounts = contra_sois_descending_pcounts[::-1]

ipsi_pcounts_descending_order = np.sort(_pcipsi_hem)
ipsi_order = np.argsort(np.mean(_pcipsi_hem, axis = 0))
ipsi_sois_descending_pcounts = np.array(short_nuclei)[ipsi_order]
ipsi_sois_ascending_pcounts = ipsi_sois_descending_pcounts[::-1]

contra_ranks = [i for i,nuc in enumerate(contra_sois_ascending_pcounts)]
ipsi_ranks = [i for j,nu in enumerate(ipsi_sois_ascending_pcounts) 
    for i,nuc in enumerate(contra_sois_ascending_pcounts) if ipsi_sois_ascending_pcounts[j] == nuc]

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#size of scatter
size = 70

Y = contra_ranks
X = ipsi_ranks

results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
ax.scatter(y = Y, x = X, s = size)

#plot fit line
ax.plot(mean_slope*range(len(X))+mean_intercept, '--k')    
    
ax.set_xlabel("Contralateral thalamus rank order (hemisphere injections)")
ax.set_ylabel("Ipsilateral thalamus rank order (hemisphere injections)")

#make text box
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)


plt.savefig(os.path.join(fig_dst, "thal_contra_v_ipsi_rank_order_hemisphere_threshold_%s.pdf" % threshold), bbox_inches = "tight")

plt.close()

#%%
#regression of total thalamic density over contra/ipsi ratio

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#size of scatter
size = 70

_Y = _ratio[1]
#filter out zeros
Y = _Y[_Y > 0]
X = (density_per_brain_left[:,0]+density_per_brain_right[:,0])[_Y > 0]

results = sm.OLS(np.log10(Y),sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
ax.scatter(y = np.log10(Y), x = X, s = size)

#plot fit line
#ax.plot(mean_slope*range(len(X))+mean_intercept, '--k')    
    
ax.set_xlabel("SM thalamic density (Cells/$mm^3$)")
ax.set_ylabel("Contra/Ipsi ratio (log 10)")
#make text box
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

textstr = "\n".join((
    "slope: {:0.5f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)

plt.savefig(os.path.join(fig_dst, "sm_thal_contra_ipsi_ratio_v_density.png"), bbox_inches = "tight")

plt.close()

#%%

#large panel figures 

## display
fig, axes = plt.subplots(ncols = 1, nrows = 5, figsize = (10,4), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,0.8,0.8,0.8,0.5]})


#set colormap specs
vmaxcounts = 70
whitetext = 10
yaxis_label_x,yaxis_label_y = -0.3,0.5

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
show = sort_dcontra.T
yaxis = ["Pons", "Sensory-motor thalamus", "Polymodal association thalamus"]

vmin = 0
vmax = vmaxcounts
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
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra", fontsize="small")
ax.yaxis.set_label_coords(yaxis_label_x,yaxis_label_y)

ax = axes[2]
show = sort_dipsi.T

vmin = 0
vmax = vmaxcounts
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
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Ipsi", fontsize="small")
ax.yaxis.set_label_coords(yaxis_label_x,yaxis_label_y)


ax = axes[3]
show = sort_ratio.T

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
cb.set_label("Ratio", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col > 1.5:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")

# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")#plt.savefig(os.path.join(dst, "thal_glm.pdf"), bbox_inches = "tight")
ax.set_ylabel("Contra/Ipsi", fontsize="small")
ax.yaxis.set_label_coords(yaxis_label_x,yaxis_label_y)

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
ax.set_xticks(np.arange(len(sort_brains))+.5)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=5, ha="right")

ax.set_yticks(np.arange(1)+.5)
ax.set_yticklabels(["M-L distance"], fontsize="x-small")

plt.savefig(os.path.join(fig_dst, "thal_contra_ipsi_ratio_w_density.pdf"), bbox_inches = "tight")

#%%
## display
fig, axes = plt.subplots(ncols = 1, nrows = 5, figsize = (10,4), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,0.8,0.8,0.8,0.5]})


#set colormap specs
vmaxcounts = 800
whitetext = 100
yaxis_label_x,yaxis_label_y = -0.3,0.5 #label coords position
    
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
show = sort_ccontra.T

vmin = 0
vmax = vmaxcounts
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
        if col < whitetext:
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
ax.yaxis.set_label_coords(yaxis_label_x,yaxis_label_y)

ax = axes[2]
show = sort_cipsi.T

vmin = 0
vmax = vmaxcounts
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
        if col < whitetext:
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
ax.yaxis.set_label_coords(yaxis_label_x,yaxis_label_y)


ax = axes[3]
show = sort_ratio.T

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
cb.set_label("Ratio", fontsize="x-small", labelpad=3)
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
ax.yaxis.set_label_coords(yaxis_label_x,yaxis_label_y)

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

plt.savefig(os.path.join(fig_dst, "thal_contra_ipsi_ratio_w_counts.pdf"), bbox_inches = "tight")

#%%
#basic statistics for these ratios
from scipy.stats import median_absolute_deviation as mad

df = pd.DataFrame()
#decimal to round by
d = 2
df["median"] = np.round(np.median(sort_ratio, axis = 0), d)
df["mean"] = np.round(np.mean(sort_ratio, axis = 0), d)
df["std"] = np.round(np.std(sort_ratio, axis = 0), d)
df["est std"] = np.round(mad(sort_ratio, axis = 0)/0.6745, d)

df.index = yaxis

df.to_csv(os.path.join(fig_dst, "thal_contra_ipsi_ratio_stats.csv"))

#%%
#regress pons ratios against contra/ipsi rations for sensory-motor thalamus

import statsmodels.api as sm

pons_ratio = _ratio[0,:]
sm_thal_ratio = _ratio[1,:]
poly_thal_ratio = _ratio[2,:]

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([.4,.1,.5,.8])

#size of scatter
size = 70

Y = poly_thal_ratio
X = pons_ratio

results = sm.OLS(Y,sm.add_constant(X)).fit()

mean_slope = results.params[1]
mean_r2 = results.rsquared
mean_intercept = results.params[0]

#plot as scatter   
ax.scatter(y = Y, x = X, s = size)

#plot fit line
#ax.plot(mean_slope*range(5)+mean_intercept, '--k')    
    
ax.set_ylim([0, 4])
ax.set_xlabel("Contra/Ipsi ratio (Pons)")
ax.set_ylabel("Contra/Ipsi ratio (Polymodal association thalamus)")

#make text box
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

ax.legend(loc="upper left")
textstr = "\n".join((
    "slope: {:0.2f}".format(mean_slope),
    "$R^2$: {:0.2f}".format(mean_r2)))

ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)


plt.savefig(os.path.join(fig_dst, "pons_v_sm_thal.pdf"), bbox_inches = "tight")