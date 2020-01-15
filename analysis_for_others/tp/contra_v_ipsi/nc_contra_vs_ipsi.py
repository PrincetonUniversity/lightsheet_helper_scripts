#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 09:17:20 2019

@author: wanglab
"""

%matplotlib inline
import numpy as np, pandas as pd, os, sys, shutil, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl, json, statsmodels.api as sm
from tools.registration.register import elastix_command_line_call, jacobian_command_line_call, change_interpolation_order, transformix_command_line_call, transformed_pnts_to_allen_helper_func, count_structure_lister
from tools.utils.io import listdirfull, makedir, load_memmap_arr, load_np, listall, load_kwargs
from skimage.external import tifffile
from tools.analysis.network_analysis import make_structure_objects
from scipy.ndimage.measurements import center_of_mass

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#set paths
dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/"
fig_dst = "/home/wanglab/Desktop"

#USING 60um edge erosion and 80 um ventricular erosion for NC, as edge seems to be the break bpoint. No real effect for ventricular so will keep the same
ann_pth = "/jukebox/wang/pisano/Python/atlas/stepwise_erosion/annotation_sagittal_atlas_20um_iso_60um_edge_erosion_80um_ventricular_erosion.tif"#"/jukebox/wang/pisano/Python/atlas/annotation_sagittal_atlas_20um_iso.tif"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx"

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
   
#collect 
id_table = pd.read_excel(df_pth)
#brains should be in this order as they were saved in this order for inj analysis
brains = ["20180409_jg46_bl6_lob6a_04","20180608_jg75","20170204_tp_bl6_cri_1750r_03","20180608_jg72",
 "20180416_jg56_bl6_lob8_04","20170116_tp_bl6_lob45_ml_11","20180417_jg60_bl6_cri_04","20180410_jg52_bl6_lob7_05",
 "20170116_tp_bl6_lob7_1000r_10","20180409_jg44_bl6_lob6a_02","20180410_jg49_bl6_lob45_02","20180410_jg48_bl6_lob6a_01",
 "20180612_jg80","20180608_jg71","20170212_tp_bl6_crii_1000r_02","20170115_tp_bl6_lob6a_rpv_03","20170212_tp_bl6_crii_2000r_03",
 "20180417_jg58_bl6_sim_02","20170130_tp_bl6_sim_1750r_03","20170115_tp_bl6_lob6b_ml_04","20180410_jg50_bl6_lob6b_03",
 "20170115_tp_bl6_lob6a_1000r_02","20170116_tp_bl6_lob45_500r_12","20180612_jg77","20180612_jg76",
 "20180416_jg55_bl6_lob8_03","20170115_tp_bl6_lob6a_500r_01","20170130_tp_bl6_sim_rpv_01","20170204_tp_bl6_cri_1000r_02",
 "20170212_tp_bl6_crii_250r_01","20180417_jg61_bl6_crii_05","20170116_tp_bl6_lob7_ml_08","20180409_jg47_bl6_lob6a_05"]
    
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/neocortex/injection_sites"
imgs = [os.path.join(src, xx+".tif.tif") for xx in brains]

#pool brain names and L/R designation into dict
lr_dist = {}
inj_vox = {}

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
    lr_dist[brain[:-4]] = dist
    inj_vox[brain[:-4]] = inj_vol
    
    if dist < 0:
        print("brain {} has a left-sided injection\n".format(brain))
    elif dist > 0:
        print("brain {} has a right-sided injection\n".format(brain))
    else:
        print("brain has an injection close to midline so not considering it rn\n")


#get injection fractions
inj_raw = np.array([inj.astype(bool) for inj in inj_vox.values()])
    
atl_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif")[:, 423:, :] #cropping coz tom segmented this
ann_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif")[:, 423:, :]
anns = np.unique(ann_raw).astype(int)
print(ann_raw.shape)
     
#annotation IDs of the cerebellum ONLY that are actually represented in annotation file
iids = {"Lingula (I)": 912,
        "Lobule II": 976,
        "Lobule III": 984,
        "Lobule IV-V": 1091,
        "Lobule VIa": 936,
        "Lobule VIb": 1134,
        "Lobule VII": 944,
        "Lobule VIII": 951,
        "Lobule IX": 957, #uvula IX
        "Lobule X": 968, #nodulus X
        "Simplex lobule": 1007, #simplex
        "Crus 1": 1056, #crus 1
        "Crus 2": 1064, #crus 2
        "Paramedian lobule": 1025, #paramedian lob
        "Copula pyramidis": 1033 #copula pyramidis
        }
ak = np.asarray([k for k,v in iids.items()])

atlas_rois = {}
for nm, iid in iids.items():
    z,y,x = np.where(ann_raw == iid) #find where structure is represented
    ann_blank = np.zeros_like(ann_raw)
    ann_blank[z,y,x] = 1 #make a mask of structure in annotation space
    atlas_rois[nm] = ann_blank.astype(bool) #add mask of structure to dictionary

expr_all_as_frac_of_lob = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
expr_all_as_frac_of_inj = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / mouse.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
primary = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj])
primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])
secondary = np.array([np.argsort(e)[-2] for e in expr_all_as_frac_of_inj])

#pooled injections
ak_pool = np.array(["Lob. I-III, IV-V", "Lob. VIa, VIb, VII", "Lob. VIII, IX, X", #no simpplex injections
                 "Simplex", "Crus I", "Crus II", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
#get n's after pooling
primary_lob_n = np.array([len(np.where(primary_pool == i)[0]) for i in range(max(primary_pool)+1)])
#normalization  of inj site
frac_of_inj_pool_norm = np.asarray([brain/brain.sum() for brain in frac_of_inj_pool])

#%%

#get brains that we actually need to get cell counts from
lr_brains = list(lr_designation.keys())
src = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts/transformed_points"

post_transformed = [os.path.join(src, os.path.join(xx, "transformed_points")) for xx in lr_brains]

#def transformed_cells_to_ann(fld, ann, dst, fl_nm):
#    """ consolidating to one function bc then no need to copy/paste """
#    dct = {}
#    
#    for fl in fld:
#        converted_points = os.path.join(fl, "posttransformed_zyx_voxels.npy")
#        print(converted_points)
#        point_lst = transformed_pnts_to_allen_helper_func(np.load(converted_points), ann, order = "ZYX")
#        df = count_structure_lister(id_table, *point_lst).fillna(0)
#        #for some reason duplicating columns, so use this
#        nm_cnt = pd.Series(df.cell_count.values, df.name.values).to_dict()
#        fl_name = os.path.basename(os.path.dirname(fl))
#        dct[fl_name]= nm_cnt
#        
#    #unpack
#    index = dct[list(dct.keys())[0]].keys()
#    columns = dct.keys()
#    data = np.asarray([[dct[col][idx] for idx in index] for col in columns])
#    df = pd.DataFrame(data.T, columns=columns, index=index)
#    
#    #save before adding projeny counts at each level
#    df.to_pickle(os.path.join(dst, fl_nm))
#    
#    return os.path.join(dst, fl_nm)
#
##collect counts from right side
#right = transformed_cells_to_ann(post_transformed, ann_right, dst, "nc_right_side_no_prog_at_each_level_pma.p")
##collect counts from left side
#left = transformed_cells_to_ann(post_transformed, ann_left, dst, "nc_left_side_no_prog_at_each_level_pma.p")


#import dict of cells by region
r_cells_regions = pckl.load(open(os.path.join(dst, "nc_right_side_no_prog_at_each_level_pma.p"), "rb"), encoding = "latin1")
r_cells_regions = r_cells_regions.to_dict(orient = "dict")      

contra = {}; ipsi = {} #collect contra and ipsi frame
for k,v in r_cells_regions.items():
    if lr_dist[k] < 0:
        contra[k] = v
    else:
        ipsi[k] = v

#LEFT SIDE
l_cells_regions = pckl.load(open(os.path.join(dst, "nc_left_side_no_prog_at_each_level_pma.p"), "rb"), encoding = "latin1")
l_cells_regions = l_cells_regions.to_dict(orient = "dict")      

for k,v in l_cells_regions.items():
    if lr_dist[k] > 0:
        contra[k] = v
    else:
        ipsi[k] = v
        
contra_df = pd.DataFrame(contra)
contra_df.to_csv(os.path.join(dst, "data/nc_contra_counts_33_brains_pma.csv")) 

ipsi_df = pd.DataFrame(ipsi)
ipsi_df.to_csv(os.path.join(dst, "data/nc_ipsi_counts_33_brains_pma.csv"))         
#%%
cells_regions_pth = os.path.join(dst, "data/nc_contra_counts_33_brains_pma.csv")

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.020
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

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
    if 'msg' in list(dic.keys()): dic = dic['msg'][0]
    
    name = dic.get('name')
    children = dic.get('children')
    if name == parent_structure:
        for child in children: # child is a dict
            child_name = child.get('name')
            progeny_list.append(child_name)
            get_progeny(child,parent_structure=child_name,progeny_list=progeny_list)
        return
    
    for child in children:
        child_name = child.get('name')
        get_progeny(child,parent_structure=parent_structure,progeny_list=progeny_list)
    return 

#get progeny of all large structures
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

#get counts for all of neocortex
sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

#%%
pcounts = np.array([xx/sum(xx) for xx in counts_per_struct.T])*100

#make % counts map like the h129 dataset (nc only for now)

## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (5,6), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [2,5]})

#set colorbar features 
maxpcount = 40
whitetext = 3
annotation_size = "x-small" #for the number annotations inside the heatmap
brain_lbl_size = "x-small"
yaxis = sois #for density by nc areas map

#sort inj fractions by primary lob
sort_pcounts = [pcounts[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_pcounts = np.array(list(itertools.chain.from_iterable(sort_pcounts)))
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

vmin = 0
vmax = 0.8
cmap = plt.cm.Reds 
cmap.set_over('darkred')
#colormap
bounds = np.linspace(vmin,vmax,4)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.7, aspect=8)
cb.set_label("% of injection", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")

ax = axes[1]
show = np.fliplr(sort_pcounts).T

vmin = 0
vmax = maxpcount
cmap = plt.cm.viridis
cmap.set_over("orange")
#colormap
bounds = np.linspace(vmin,vmax,((vmax-vmin)/5)+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, 
                  format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("% of total neocortical counts", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(np.flipud(yaxis), fontsize="x-small")

ax.set_xticks(np.arange(len(sort_brains))+.5)
lbls = np.asarray(sort_brains)
ax.set_xticklabels(sort_brains, rotation=30, fontsize=brain_lbl_size, ha="right")

plt.savefig(os.path.join(fig_dst, "h129_pcounts_nc.pdf"), bbox_inches = "tight")

#%%
pcounts_pool = np.asarray([[brain[0]+brain[1]+brain[2]+brain[4], brain[3], brain[6], brain[5]+brain[7], 
                                          brain[8]+brain[9], brain[10], brain[12], brain[11], 
                                          brain[13]+brain[14], brain[15]+brain[16]] for brain in pcounts])


regions = np.asarray(["Infralimbic, Prelimbic,\n Ant. Cingulate, Orbital",
       "Frontal pole", "Agranular insula", "Gustatory, Visceral",
       "Somatomotor, Somatosensory", "Retrosplenial", "Visual",
       "Post. parietal", "Temporal, Auditory", "Peririhinal, Ectorhinal"])
    
    
X = frac_of_inj_pool_norm
Y = pcounts_pool    

#%%

##  glm
c_mat = []
mat = []
pmat = []
mat_shuf = []
p_shuf = []
fit = []
fit_shuf = []

for itera in range(1000):
    if itera%100 == 0: print(itera)
    if itera == 0:
        shuffle = False
        inj = X.copy()
    else:
        shuffle = True
        mat_shuf.append([])
        p_shuf.append([])
        fit_shuf.append([])
        inj = X[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:]
    for count, region in zip(Y.T, regions):
        inj_ = inj[~np.isnan(count)]
        count = count[~np.isnan(count)]

        # intercept:
        inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
        
#        glm = sm.OLS(count, inj_)
        glm = sm.GLM(count, inj_, family=sm.families.Poisson())
        res = glm.fit()
        
        coef = res.params[:-1]
        se = res.bse[:-1]
        pvals = res.pvalues[:-1] 

        val = coef/se

        if not shuffle:
            c_mat.append(coef)
            mat.append(val)
            pmat.append(pvals)
            fit.append(res.fittedvalues)
            
        elif shuffle:
            mat_shuf[-1].append(val)
            p_shuf[-1].append(pvals)
            fit_shuf[-1].append(res.fittedvalues)
        # inspect residuals
        """
        if not shuffle:
            plt.clf()
            plt.scatter(res.fittedvalues, res.resid)
            plt.hlines(0, res.fittedvalues.min(), res.fittedvalues.max())
            plt.title(region)
            plt.xlabel("Fit-vals")
            plt.ylabel("Residuals")
            plt.savefig(os.path.join(dst, "resid_inspection-{}.png").format(region))
        """
        
c_mat = np.array(c_mat)
mat = np.array(mat) # region x inj
mat_shuf = np.array(mat_shuf) # region x inj
pmat = np.array(pmat) # region x inj
p_shuf = np.array(p_shuf)
fit = np.array(fit)
fit_shuf = np.array(fit_shuf)

#%%

## display
fig = plt.figure(figsize=(8,5))
ax = fig.add_axes([.4,.1,.5,.8])

# map 1: weights
show = np.flipud(mat) # NOTE abs

vmin = 1
vmax = 6
whitetext = 4
fontsize = "medium"
cmap = plt.cm.Reds
cmap.set_under("w")
cmap.set_over("maroon")
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,(vmax-vmin) + 1)
#bounds = np.linspace(0,5,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", ticks=bounds, boundaries=bounds, format="%0.1f", 
                  shrink=0.3, aspect=10)
cb.set_label("Weight / SE", fontsize="x-small", labelpad=3)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)

# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        if col < whitetext:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize=fontsize)
        else: 
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize=fontsize)

# signif
sig = np.flipud(pmat) < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*10)
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x, y+0.3, "*", fontsize=10, ha="left", va="bottom", color = "black", transform=ax.transData)
#ax.text(.5, 1.06, "X: p<0.05".format(nullmean, nullstd), ha="center", va="center", fontsize="small", transform=ax.transAxes)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
# xticksjt -t monokai -m 200
ax.set_xticks(np.arange(len(ak_pool))+.5)

#remaking labeles so it doesn"t look squished
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], rotation=30, fontsize="x-small", ha="right")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(["{}".format(bi) for bi in np.flipud(regions)], fontsize=fontsize)
plt.savefig(os.path.join(dst, "h129_nc_glm.pdf"), bbox_inches = "tight")