#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 18:57:00 2019

@author: wanglab
"""


import matplotlib as mpl, os, pandas as pd, json, statsmodels.api as sm
import matplotlib.pyplot as plt
import numpy as np, pickle as pckl

#custom
src = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
cells_regions_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/nc_contra_counts_33_brains_pma.csv"
dst = "/home/wanglab/Desktop/"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

data_pth = os.path.join(src, "nc_hsv_maps_contra_pma.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

#set the appropritate variables
brains = data["brains"]
expr_all_as_frac_of_inj = data["expr_all_as_frac_of_inj"]
ak_pool = data["ak_pool"]
frac_of_lob = data["expr_all_as_frac_of_lob"]

brains = ['20180409_jg46_bl6_lob6a_04', '20180608_jg75',
       '20170204_tp_bl6_cri_1750r_03', '20180608_jg72',
       '20180416_jg56_bl6_lob8_04', '20170116_tp_bl6_lob45_ml_11',
       '20180417_jg60_bl6_cri_04', '20180410_jg52_bl6_lob7_05',
       '20170116_tp_bl6_lob7_1000r_10', '20180409_jg44_bl6_lob6a_02',
       '20180410_jg49_bl6_lob45_02', '20180410_jg48_bl6_lob6a_01',
       '20180612_jg80', '20180608_jg71', '20170212_tp_bl6_crii_1000r_02',
       '20170115_tp_bl6_lob6a_rpv_03', '20170212_tp_bl6_crii_2000r_03',
       '20180417_jg58_bl6_sim_02', '20170130_tp_bl6_sim_1750r_03',
       '20170115_tp_bl6_lob6b_ml_04', '20180410_jg50_bl6_lob6b_03',
       '20170115_tp_bl6_lob6a_1000r_02', '20170116_tp_bl6_lob45_500r_12',
       '20180612_jg77', '20180612_jg76', '20180416_jg55_bl6_lob8_03',
       '20170115_tp_bl6_lob6a_500r_01', '20170130_tp_bl6_sim_rpv_01',
       '20170204_tp_bl6_cri_1000r_02', '20170212_tp_bl6_crii_250r_01',
       '20180417_jg61_bl6_crii_05', '20170116_tp_bl6_lob7_ml_08',
       '20180409_jg47_bl6_lob6a_05']
#%%

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.025
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
    #add counts from main structure
    counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)      

pcounts = np.nan_to_num(np.asarray([((brain/sum(brain))*100) for brain in counts_per_struct.T]))    

density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T
#%%

pcounts = np.nan_to_num(np.asarray([((brain/sum(brain))*100) for brain in counts_per_struct.T]))    

pcounts_pool = np.asarray([[xx[0]+xx[1]+xx[2]+xx[4], xx[3], xx[6], xx[5]+xx[7], 
                                          xx[8]+xx[9], xx[10], xx[12], xx[11], 
                                          xx[13]+xx[14], xx[15]+xx[16]] for xx in pcounts])

vol_pool = np.asarray([vol[0]+vol[1]+vol[2]+vol[4], vol[3], vol[6], vol[5]+vol[7], 
                                          vol[8]+vol[9], vol[10], vol[12], vol[11], 
                                          vol[13]+vol[14], vol[15]+vol[16]])

#can fix this later
sois_pool = np.asarray([sois[0]+sois[1]+sois[2]+sois[4], sois[3], sois[6], sois[5]+sois[7], 
                                          sois[8]+sois[9], sois[10], sois[12], sois[11], 
                                          sois[13]+sois[14], sois[15]+sois[16]])
#sort pcount groups by region size
sois_pool = sois_pool[np.argsort(vol_pool)]
pcounts_pool = pcounts_pool.T[np.argsort(vol_pool)].T

regions = np.array(['F Pole',
                   'P Par',
                   'Pr, EcR', 'Gust, Visc',
                   'Insula',
                   'Temp, Aud', 'RS',
                   'VIS',
                   'IL, PrL,\nAC, Orb',
                   'SM, SS'])
    
#pooled injections
ak_pool = np.array(["Lob. I-III, IV-V", "Lob. VIa, VIb, VII", "Lob. VIII, IX, X", #no simpplex injections
                 "Simplex", "Crus I", "Crus II", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
frac_of_inj_pool_norm = np.array([xx/sum(xx) for xx in frac_of_inj_pool])    

X = frac_of_inj_pool_norm
Y = pcounts_pool    

##try without pooled NC regions
pcounts = pcounts.T[np.argsort(vol)].T
regions = np.array(["IL", "PrL", "AC", "F Pole", "Orb", "Gust", "Insula", "Visc", "SM", "SS", "RS", "P Par", "VIS", 
                    "Temp", "Aud", "EcR", "Pr"]) 
regions = regions[np.argsort(vol)]

Y = pcounts

##try with density?
frac_of_lob_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in frac_of_lob])

#sort density
density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T
density = (density.T[np.argsort(vol)]).T
Y = density

#normalize by density of cells in all of isocortex?
#get counts for all of neocortex
reg = ["Isocortex"]

#first calculate counts across entire nc region
iscortex = []
for soi in reg:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    iscortex.append(np.array(counts).sum(axis = 0))
iscortex = np.array(iscortex)
iscortex_vol = []
for soi in reg:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    for progen in progeny:
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    iscortex_vol.append(np.array(counts).sum(axis = 0))
iscortex_vol = np.array(vol)  

isocortex_density = np.array([xx/(iscortex_vol[i]*(scale_factor**3)) for i, xx in enumerate(iscortex)]).T

density_norm = (density/isocortex_density)

X = frac_of_lob_pool

Y = density

#%%
#check out mean density per injection site?

#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(density[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([.4,.1,.5,.8])

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 200
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))
#colormap

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.3, aspect=10)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(lbls, rotation="vertical")
# yticks
ax.set_yticks(np.arange(len(regions))+.5)

ax.set_yticklabels(regions, fontsize="small")

#%%

##  glm
c_mat = []
mat = []
pmat = []
mat_shuf = []
p_shuf = []
fit = []
fit_shuf = []

for itera in range(100):
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
fig,ax = plt.subplots(figsize=(3,6))

# map 1: weights
show = mat

# SET COLORMAP
vmin = 0
vmax = 15
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))
annotation = False

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("Model weight / SE", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# exact value annotations
if annotation:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            if col > 5:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="small")
            else:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="small")

# signif
#only get positive significance??             
pmat_pos = np.where(mat > 0, pmat, pmat*np.inf)
sig = pmat_pos < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*np.inf)
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x+0.5, y+0.4, "*", fontsize=12, horizontalalignment='center', verticalalignment='center',
            color = "k", transform=ax.transData)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)
ax.set_xticklabels(ak_pool, rotation="vertical", fontsize=10)

# yticks
ax.set_yticks(np.arange(len(regions))+.5)
ax.set_yticklabels(regions, fontsize=10)
ax.tick_params(length=6)

plt.savefig(os.path.join(dst, "hsv_nc_density_glm_contra_pma.pdf"), bbox_inches = "tight")

plt.close()