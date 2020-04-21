#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 18:57:00 2019

@author: wanglab
"""

import matplotlib as mpl, os, pandas as pd, json, statsmodels.api as sm, statsmodels.formula.api as smf
from patsy import dmatrices
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
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6

data_pth = os.path.join(src, "nc_hsv_maps_contra_pma.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

mask = [True]*33#to hide high count brain
mask[12] = False

#set the appropritate variables
expr_all_as_frac_of_inj = data["expr_all_as_frac_of_inj"][mask]
ak_pool = data["ak_pool"]
frac_of_lob = data["expr_all_as_frac_of_lob"][mask]

brains = np.array(data["brains"])[mask]
#%%

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
sois = ["Somatosensory areas", "Somatomotor areas", "Visual areas",
       "Retrosplenial area", "Agranular insular area", "Auditory areas",
       "Anterior cingulate area", "Orbital area",
       "Temporal association areas",
       "Posterior parietal association areas", "Prelimbic area",
       "Visceral area", "Ectorhinal area", "Gustatory areas",
       "Perirhinal area", "Infralimbic area",
       "Frontal pole, cerebral cortex"]

#first calculate counts across entire nc region
counts_per_struct = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    except:
        print(soi)
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

vol = []
for soi in sois:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0]/2)
    except:
        print(soi)
    for progen in progeny:
            counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0]/2)
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)      

pcounts = np.nan_to_num(np.asarray([((brain/sum(brain))*100) for brain in counts_per_struct.T]))    

density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T
#%%

pcounts = np.nan_to_num(np.asarray([((brain/sum(brain))*100) for brain in counts_per_struct.T]))    
    
#pooled injections
ak_pool = np.array(["Lob. I-III, IV-V", "Lob. VIa, VIb, VII", "Lob. VIII, IX, X", #no simpplex injections
                 "Simplex", "Crus I", "Crus II", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
frac_of_inj_pool_norm = np.array([xx/sum(xx) for xx in frac_of_inj_pool])    

##try with density?
frac_of_lob_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]), xx[10], xx[11], xx[12], np.sum(xx[13:16])] 
                                for xx in frac_of_lob])
primary_lob_n = np.array([len(np.where(primary_pool == i)[0]) for i in range(max(primary_pool)+1)])

primary_lob_pool = np.array([np.argmax(e) for e in frac_of_lob_pool])
density = np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T

#%%

#check out mean density per injection site?
#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.median(density[np.where(primary_lob_pool == idx)[0]], axis=0) for idx in np.unique(primary_lob_pool)])

fig, ax = plt.subplots(figsize=(3.3,6))

show = mean_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 400
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))
#colormap

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.3, aspect=10)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(lbls, rotation="vertical")
# yticks
ax.set_yticks(np.arange(len(sois))+.5)

ax.set_yticklabels(sois, fontsize="small")

plt.savefig(os.path.join(dst, "hsv_nc_mean_density_contra_pma.pdf"), bbox_inches = "tight")

#%%

# X = np.array([(d-d.min())/(d.max()-d.min()) for d in frac_of_lob_pool]) #normalized... frac_of_lob_pool
# Y = np.array([(d-d.min())/(d.max()-d.min()) for d in density]) #normalized...

X = frac_of_inj_pool_norm
Y = density
#%%
#fit poisson model to find mu
#https://towardsdatascience.com/negative-binomial-regression-f99031bb25b4

alpha = []
for i,j in enumerate(Y.T):
    
    print(sois[i])
    df = pd.DataFrame(X)
    df.columns = ["l1_3", "l6_7", "l8_10", "sim", "cr1", "cr2", "pm_cp"]
    df["Density"] = j
    
    mask = np.random.rand(len(df)) < 1
    df_train = df[mask]
    df_test = df[~mask]
    print('\nTraining data set length='+str(len(df_train)))
    print('\nTesting data set length='+str(len(df_test)))
    
    expr = """Density ~ l1_3 + l6_7 + l8_10 + sim + cr1 + cr2 + pm_cp"""
    
    y_train, X_train = dmatrices(expr, df_train, return_type='dataframe')
    y_test, X_test = dmatrices(expr, df_test, return_type='dataframe')
    
    #fit
    poisson_training_results = sm.GLM(y_train, X_train, family=sm.families.Poisson()).fit()
    
    #fit OLS to find alpha
    df_train["LAMBDA"] = poisson_training_results.mu
    
    df_train['AUX_OLS_DEP'] = df_train.apply(lambda x: ((x['Density'] - x['LAMBDA'])**2 - x['Density']) / x['LAMBDA'], 
                                             axis=1)
    
    ols_expr = """AUX_OLS_DEP ~ LAMBDA - 1"""
    
    #fit ols
    aux_olsr_results = smf.ols(ols_expr, df_train).fit()
    print(aux_olsr_results.params.LAMBDA, aux_olsr_results.pvalues)
    
    alpha.append(aux_olsr_results.params.LAMBDA)
    
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
        
    i = 0
    for count, region in zip(Y.T, sois):
        inj_ = inj[~np.isnan(count)]
        count = count[~np.isnan(count)]

        # intercept:
        inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
        
        # glm = sm.OLS(count, inj_)
        glm = sm.GLM(count, inj_, family=sm.families.NegativeBinomial(alpha = alpha[i])) 
        res = glm.fit()
        
        coef = res.params[:-1]
        se = res.bse[:-1]
        pvals = res.pvalues[:-1] 

        val = coef/se

        if not shuffle:
            print(region, val)
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
        i += 1
        
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
show = np.flipud(mat) #can we use abs???

# SET COLORMAP
vmin = 0
vmax = 3.5
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
sig = np.flipud(pmat) < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*np.inf)
null = (p_shuf < .05).sum(axis=(1,2))
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
ax.set_yticks(np.arange(len(sois))+.5)
ax.set_yticklabels(np.flipud(sois), fontsize=10)

plt.savefig(os.path.join(dst, "hsv_nc_density_glm_contra_pma.pdf"), bbox_inches = "tight")
plt.close()