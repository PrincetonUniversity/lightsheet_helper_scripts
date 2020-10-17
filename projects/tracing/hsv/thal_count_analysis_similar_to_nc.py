#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 17:24:45 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl
import matplotlib as mpl, json, itertools, statsmodels.api as sm, copy
from patsy import dmatrices
import matplotlib.colors, statsmodels.formula.api as smf
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #lime color makes cells pop
 
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 3
mpl.rcParams["ytick.major.size"] = 3

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/"
fig_dst = "/home/wanglab/Desktop"

ann_pth = os.path.join(dst, "atlases/sagittal_allen_ann_25um_iso_60um_edge_80um_ventricular_erosion.tif")
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
ontology_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen.json"

#collect 
data_pth = os.path.join(dst, "data/thal_hsv_maps_contra_allen.p")
data = pckl.load(open(data_pth, "rb"), encoding = "latin1")

primary_pool = data["primary_pool"]
frac_of_inj_pool = data["frac_of_inj_pool"]
ak_pool = data["ak_pool"]
brains = np.array(data["brains"])

primary_lob_n = np.asarray([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])
frac_of_inj_pool_norm = np.asarray([brain/brain.sum() for brain in frac_of_inj_pool])

#%%
cells_regions_pth = os.path.join(dst, "data/thal_contra_counts_23_brains_80um_ventric_erosion.csv")

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

#get progeny of all large structures
with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

nuclei = ["Thalamus", "Zona incerta", "Ventral posteromedial nucleus of the thalamus",
       "Reticular nucleus of the thalamus",
       "Mediodorsal nucleus of thalamus",
       "Posterior complex of the thalamus",
       "Lateral posterior nucleus of the thalamus",
       "Lateral dorsal nucleus of thalamus",
       "Ventral medial nucleus of the thalamus",
       "Ventral posterolateral nucleus of the thalamus",
       "Ventral anterior-lateral complex of the thalamus",
       "Medial geniculate complex",
       "Dorsal part of the lateral geniculate complex",
       "Nucleus of reuniens", "Paraventricular nucleus of the thalamus",
       "Anteroventral nucleus of thalamus", "Anteromedial nucleus",
       "Parafascicular nucleus", "Ventral part of the lateral geniculate complex",
       "Medial habenula", "Central lateral nucleus of the thalamus",
       "Lateral habenula", "Submedial nucleus of the thalamus",
       "Parataenial nucleus", "Subparafascicular nucleus",
       "Central medial nucleus of the thalamus", "Anterodorsal nucleus",
       "Paracentral nucleus"]

#first calculate counts across entire nc region
counts_per_struct = []
for soi in nuclei:
    progeny = []; counts = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    counts.append([cells_regions.loc[cells_regions.Structure == soi, brain].values[0] for brain in brains])
    for progen in progeny:
        counts.append([cells_regions.loc[cells_regions.Structure == progen, brain].values[0] for brain in brains])
    counts_per_struct.append(np.array(counts).sum(axis = 0))
counts_per_struct = np.array(counts_per_struct)

pcounts = np.nan_to_num(np.asarray([((brain[1:]/brain[0])*100) for brain in counts_per_struct.T]))    

#voxels
vol = []
for soi in nuclei:
    progeny = []; counts = []; iids = []
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

density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct)]).T) #remove thalamus

#remove thalamus from density
nuclei = nuclei[1:]
density = density[:, 1:]
counts_per_struct = counts_per_struct[1:,:]
#%%

#boxplots for counts
import seaborn as sns

#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
order = np.argsort(np.median(counts_per_struct.T, axis = 0))[::-1]
sois_sort = np.array(nuclei)[order][:10]

#boxplots of percent counts
plt.figure(figsize = (5,4))
df = pd.DataFrame(pcounts)
df.columns = nuclei
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort)
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={'facecolor':'None'}, order = sois_sort)
plt.xlabel("# Neurons")
plt.ylabel("Subnucleus")
plt.savefig(os.path.join(fig_dst, "thal_counts_boxplots.pdf"), bbox_inches = "tight")

#%%

#boxplots of density
#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
order = np.argsort(np.median(density, axis = 0))[::-1]
sois_sort = np.array(nuclei)[order][:10]

#boxplots of percent counts
plt.figure(figsize = (5,4))
df = pd.DataFrame(density)
df.columns = nuclei
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort)
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={'facecolor':'None'}, order = sois_sort)
plt.xlabel("Cells / mm$^3$")
plt.ylabel("Subnucleus")
# plt.xlim([-10, 100])
plt.tick_params(length=6)

plt.savefig(os.path.join(fig_dst, "thal_density_boxplots.pdf"), bbox_inches = "tight")
#%%

#make injection site heatmap only
fig, ax = plt.subplots(figsize = (5,2))

sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
show = np.fliplr(sort_inj).T

#colormap settings
cmap = copy.copy(plt.cm.Reds)
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)#np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
lbls = np.asarray(sort_brains)
ax.set_xticklabels(np.array([ 1,  5, 10, 15, 20]))
ax.tick_params(length=6)

plt.savefig(os.path.join(fig_dst, "hsv_inj_thal_red.pdf"), bbox_inches = "tight")

#%%
#set colorbar features 
maxpcount = 8
yaxis = np.flipud(nuclei)

#make % counts map like the h129 dataset
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (5,7), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,5]})


#sort inj fractions by primary lob
sort_pcounts = [pcounts[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_pcounts = c
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

cmap = plt.cm.Reds 
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
ax.tick_params(length=6)

ax = axes[1]
show = np.fliplr(sort_pcounts).T

# SET COLORMAP
vmin = 0
vmax = maxpcount
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.4)
cb.set_label("% of thalamic neurons", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="x-small")
ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

plt.savefig(os.path.join(fig_dst, "hsv_pcounts_thal.pdf"), bbox_inches = "tight")

#%%

#set colorbar features 
maxdensity = 40
yaxis = np.flipud(nuclei)

#make density map
## display
fig, axes = plt.subplots(ncols = 1, nrows = 2, figsize = (5,7), sharex = True, gridspec_kw = {"wspace":0, "hspace":0,
                         "height_ratios": [1.5,5]})


#sort inj fractions by primary lob
sort_density = [density[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_density = np.array(list(itertools.chain.from_iterable(sort_density)))
sort_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_inj = [frac_of_inj_pool[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sort_brains = list(itertools.chain.from_iterable(sort_brains))
sort_inj = np.array(list(itertools.chain.from_iterable(sort_inj)))

#inj fractions
ax = axes[0]
show = np.fliplr(sort_inj).T

cmap = plt.cm.Reds 
cmap.set_over(cmap(1.0))
cmap.set_under("white")
vmin = 0.05
vmax = 0.8

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.8)#
cb.set_label("Injection % coverage\n of region", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True) #TP
ax.set_yticks(np.arange(len(ak_pool))+.5)
ax.set_yticklabels(np.flipud(ak_pool), fontsize="small")
ax.tick_params(length=6)

ax = axes[1]
show = np.fliplr(sort_density).T

# SET COLORMAP
vmin = 0
vmax = maxdensity
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.4)
cb.set_label("Cells / mm$^3$", fontsize="small", labelpad=5)
cb.ax.tick_params(labelsize="small")
cb.ax.set_visible(True)

# aesthetics
# yticks
ax.set_yticks(np.arange(len(yaxis))+.5)
ax.set_yticklabels(yaxis, fontsize="small")

ax.set_xticks(np.arange(0, len(sort_brains), 5)+.5)
ax.set_xticklabels(np.arange(0, len(sort_brains), 5)+1)

plt.savefig(os.path.join(fig_dst, "hsv_density_thal.pdf"), bbox_inches = "tight")
#%%
#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(pcounts[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig,ax = plt.subplots(figsize=(3,9))

show = np.flipud(mean_counts.T) 

# SET COLORMAP
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))
annotations = False

#set min and max of colorbar
vmin = 0
vmax = 8

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("Mean % of thalamic neurons", fontsize="medium", labelpad=5)
cb.ax.tick_params(labelsize="medium")

cb.ax.set_visible(True)
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], 
                   rotation="vertical", fontsize=8)
# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(np.flipud(nuclei), fontsize="medium")
ax.tick_params(length=6)

plt.savefig(os.path.join(fig_dst,"thal_mean_count.pdf"), bbox_inches = "tight")

#%%
#only look at mean counts per "cerebellar region" (i.e. that which had the highest contribution of the injection)    
mean_counts = np.asarray([np.mean(density[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

fig,ax = plt.subplots(figsize=(3,9))

show = mean_counts.T[::-1]

# SET COLORMAP
cmap = plt.cm.Blues
cmap.set_over(cmap(1.0))
annotations = False

#set min and max of colorbar
vmin = 0
vmax = 40

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, format="%d", shrink=0.3, aspect=10)
cb.set_label("Mean neurons / mm$^3$", fontsize="medium", labelpad=5)
cb.ax.tick_params(labelsize="medium")

cb.ax.set_visible(True)
        
#remaking labeles so it doesn"t look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{} ({})".format(ak, n) for ak, n in zip(lbls, primary_lob_n)], 
                   rotation="vertical", fontsize="medium")
# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(nuclei[::-1], fontsize="medium")

plt.savefig(os.path.join(fig_dst,"thal_mean_density.pdf"), bbox_inches = "tight")

#%%
#glm
X = frac_of_inj_pool_norm
Y = pcounts
# Y = np.nan_to_num(np.array([(d/sum(d))*100 for d in density]))

c_mat = []
mat = []
pmat = []
mat_shuf = []
p_shuf = []
ars = []
rs = []
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
    for count, region in zip(Y.T, nuclei):
        try:
            inj_ = inj[~np.isnan(count)]
            count = count[~np.isnan(count)]
    
            # intercept:
            inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
            
#            glm = sm.OLS(count, inj_)
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
        except:
            inj = X[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:] #sometimes the shuffle messes stuff up
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
                mat.append(val)
                pmat.append(pvals)
            elif shuffle:
                mat_shuf[-1].append(val)
                p_shuf[-1].append(pvals)
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
    

mat = np.array(mat) # region x inj
mat_shuf = np.array(mat_shuf) # region x inj
pmat = np.array(pmat) # region x inj
p_shuf = np.array(p_shuf)
fit = np.array(fit)
fit_shuf = np.array(fit_shuf)

#%%
## display
fig,ax = plt.subplots(figsize=(3,9))

# map 1: weights
show = np.flipud(mat)

vmin = 0
vmax = 4
whitetext = 4
annotation = False

# SET COLORMAP
cmap = copy.copy(plt.cm.Blues)
cmap.set_over(cmap(1.0))
annotations = False

#colormap
pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)
cb = plt.colorbar(pc, ax=ax, format="%0.1f", shrink=0.3, aspect=10)
cb.set_label("Model weight / SE", fontsize="medium", labelpad=5)
cb.ax.tick_params(labelsize="medium")
cb.ax.set_visible(True)

# exact value annotations
if annotation:
    for ri,row in enumerate(show):
        for ci,col in enumerate(row):
            if col > whitetext:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="small")
            else:
                ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="small")

# signif
pmat_pos = np.where(mat > 0, pmat, pmat*100000000000)
sig = np.flipud(pmat_pos) < .05#/np.size(pmat)
p_shuf_pos = np.where(mat_shuf < 0, p_shuf, p_shuf*10)
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x+0.5, y+0.4, "*", fontsize=18, horizontalalignment='center', verticalalignment='center',
            color = "w", transform=ax.transData)
# ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)
ax.set_xticklabels(ak_pool, rotation="vertical", fontsize="medium")

# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(np.flipud(nuclei), fontsize="medium")

plt.savefig(os.path.join(fig_dst, "thal_pcount_glm_contra_allen.pdf"), bbox_inches = "tight")

#%%
###DOES NOT WORK###
#hide high count brain for model?
mask = [True]*23
mask[6] = False
# mask[21] = False

X = frac_of_inj_pool_norm[mask]
Y = density[mask]
#fit poisson model to find mu
#https://towardsdatascience.com/negative-binomial-regression-f99031bb25b4

alpha = []
for i,j in enumerate(Y.T):
    
    print(nuclei[i])
    df = pd.DataFrame(X)
    df.columns = ["l1_3", "l6_7", "l8_10", "sim", "cr1", "cr2", "pm_cp"]
    df["Density"] = j
    
    msk = np.random.rand(len(df)) < 1
    df_train = df[msk]
    df_test = df[~msk]
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
    print(aux_olsr_results.params.LAMBDA, aux_olsr_results.pvalues.LAMBDA)
    
    alpha.append(aux_olsr_results.params.LAMBDA)
    

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
    for count, region in zip(Y.T, nuclei):
        try:
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
                c_mat.append(coef)
                mat.append(val)
                pmat.append(pvals)
                fit.append(res.fittedvalues)
                
            elif shuffle:
                mat_shuf[-1].append(val)
                p_shuf[-1].append(pvals)
                fit_shuf[-1].append(res.fittedvalues)
        except:
            inj = X[np.random.choice(np.arange(len(inj)), replace=False, size=len(inj)),:] #sometimes the shuffle messes stuff up
            inj_ = inj[~np.isnan(count)]
            count = count[~np.isnan(count)]
    
            # intercept:
            inj_ = np.concatenate([inj_, np.ones(inj_.shape[0])[:,None]*1], axis=1)
            
            glm = sm.GLM(count, inj_, family=sm.families.NegativeBinomial(alpha = alpha[i])) 
            res = glm.fit()
            
            coef = res.params[:-1]
            se = res.bse[:-1]
            pvals = res.pvalues[:-1] 
    
            val = coef/se
    
            if not shuffle:
                mat.append(val)
                pmat.append(pvals)
            elif shuffle:
                mat_shuf[-1].append(val)
                p_shuf[-1].append(pvals)
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
fig,ax = plt.subplots(figsize=(3,9))

# map 1: weights
show = np.flipud(mat) #can we use abs???

# SET COLORMAP
vmin = 0
vmax = 3
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
null = (p_shuf_pos < .05).sum(axis=(1,2))
nullmean = null.mean()
nullstd = null.std()
for y,x in np.argwhere(sig):
    pass
    ax.text(x+0.5, y+0.4, "*", fontsize=12, horizontalalignment='center', verticalalignment='center',
            color = "white", transform=ax.transData)
ax.text(.5, 1.06, "*: p<0.05\n{:0.1f} ($\pm$ {:0.1f}) *'s are expected by chance if no real effect exists".format(nullmean, nullstd), ha="center", va="center", fontsize="x-small", transform=ax.transAxes)

# aesthetics
ax.set_xticks(np.arange(len(ak_pool))+.5)
ax.set_xticklabels(ak_pool, rotation="vertical", fontsize=10)

# yticks
ax.set_yticks(np.arange(len(nuclei))+.5)
ax.set_yticklabels(np.flipud(nuclei), fontsize=10)

plt.savefig(os.path.join(fig_dst, "thal_density_glm_contra_allen.pdf"), bbox_inches = "tight")
plt.close()