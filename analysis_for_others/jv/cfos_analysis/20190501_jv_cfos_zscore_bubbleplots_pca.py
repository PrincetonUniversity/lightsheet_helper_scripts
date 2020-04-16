#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 13:50:36 2019

@author: wanglab
"""

#The Kruskal–Wallis test by ranks, Kruskal–Wallis H test[1] (named after William Kruskal and W. Allen Wallis), or one-way ANOVA on ranks[1] 
#is a non-parametric method for testing whether samples originate from the same distribution.[2][3][4] It is used for comparing two or 
#more independent samples of equal or different sample sizes. It extends the Mann–Whitney U test, which is used for comparing only two groups. 
#The parametric equivalent of the Kruskal–Wallis test is the one-way analysis of variance (ANOVA). 
import pandas as pd
from scipy.stats import kruskal
import numpy as np, os

src = "/home/wanglab/Desktop"

df_pth = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/summed_parents_cell_counts_dataframe.csv"
df = pd.read_csv(df_pth, index_col = None)
#temp change of column name
df["name"] = df["parent_name"]

df_anova = pd.DataFrame()
df_anova["name"] = np.unique(df["name"].values)

for nm in np.unique(df.name.values): #only gets unique names
    f, pval = kruskal(df[(df.name == nm) & (df.Condition == "DREADDs")].counts.values, 
                 df[(df.name == nm) & (df.Condition == "homecage_control")].counts.values,
                 df[(df.name == nm) & (df.Condition == "CNO_control_no_reversal")].counts.values, 
                 df[(df.name == nm) & (df.Condition == "CNO_control_reversal") & (df.Brain != "an17")].counts.values)
    
    df_anova.loc[(df_anova["name"] == nm), "kruskal_counts_pval"] = pval
        
    f, pval = kruskal(df[(df.name == nm) & (df.Condition == "DREADDs")]["percent"].values, 
                 df[(df.name == nm) & (df.Condition == "homecage_control")]["percent"].values,
                 df[(df.name == nm) & (df.Condition == "CNO_control_no_reversal")]["percent"].values, 
                 df[(df.name == nm) & (df.Condition == "CNO_control_reversal") & (df.Brain != "an17")]["percent"].values)

    df_anova.loc[(df_anova["name"] == nm), "kruskal_percent_counts_pval"] = pval
    
        
df_anova.to_csv(os.path.join(src, "kruskal_pooled_structures.csv"), index = None)

#%%
all_pth = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/cell_counts_dataframe_w_percents_density.csv"

df = pd.read_csv(all_pth, index_col = None)

df_anova = pd.DataFrame()
df_anova["name"] = np.unique(df["name"].values)

for nm in np.unique(df.name.values): #only gets unique names
    
    try:
        f, pval = kruskal(df[(df.name == nm) & (df.Condition == "DREADDs")].counts.values, 
                     df[(df.name == nm) & (df.Condition == "homecage_control")].counts.values,
                     df[(df.name == nm) & (df.Condition == "CNO_control_no_reversal")].counts.values, 
                     df[(df.name == nm) & (df.Condition == "CNO_control_reversal") & (df.Brain != "an17")].counts.values)
        
        df_anova.loc[(df_anova["name"] == nm), "kruskal_counts_pval"] = pval
            
        f, pval = kruskal(df[(df.name == nm) & (df.Condition == "DREADDs")]["percent"].values, 
                     df[(df.name == nm) & (df.Condition == "homecage_control")]["percent"].values,
                     df[(df.name == nm) & (df.Condition == "CNO_control_no_reversal")]["percent"].values, 
                     df[(df.name == nm) & (df.Condition == "CNO_control_reversal") & (df.Brain != "an17")]["percent"].values)
        
        df_anova.loc[(df_anova["name"] == nm), "kruskal_percent_counts_pval"] = pval
    except:
        print(nm)
        
df_anova.to_csv(os.path.join(src, "kruskal_all_structures.csv"), index = None)

#%%

import math
import matplotlib.pyplot as plt
from lightsheet.network_analysis import make_structure_objects

#bubble plot
df1 = pd.read_csv("/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/parents_w_zscores_per_structure.csv")
df2 = pd.read_csv("/home/wanglab/Desktop/kruskal_pooled_structures.csv")
df3 = pd.read_excel("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx")
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"

#which ones not to visualise
to_exclude = ["Nucleus of the lateral olfactory tract",
              "Cerebral cortex", "Cerebral nuclei", 
              "Inferior colliculus", "Perihypoglossal nuclei",
              "Cortical subplate"]
#look at stuff per condition
sig_str = df2[(df2.kruskal_percent_counts_pval < 0.05) & (~df2.name.isin(to_exclude))].name.values

#get structure ontology
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

#%%

#dreadds
x = []; s = []; abbrev = []; voxels = []
for soi in sig_str:
    x.append(np.mean(df1[(df1.parent_name == soi) & (df1.Condition == "DREADDs")]["z_score_percents"].values))
    soi = [s for s in structures if s.name==soi][0]
    print(soi.name)
    if not df3.loc[df3.name == soi.name,"voxels_in_structure"].values == 0:
        voxels.append(df3.loc[df3.name == soi.name,"voxels_in_structure"].values)
    progeny = [str(xx.name) for xx in soi.progeny]
    if len(progeny) > 0:
        for progen in progeny:
            voxels.append(df3.loc[df3.name == progen,"voxels_in_structure"].values)
    s.append(np.sum(np.asarray(voxels)))
    abbrev.append(df3[df3.name == soi.name]["acronym"].values) 

size = [xx*(0.020**3)/10 for xx in s]
acronym = [iid[0] for iid in abbrev]

y = df2[df2.name.isin(sig_str)].kruskal_percent_counts_pval.values
y = np.asarray([-math.log(xx, 10) for xx in y])

df_dreadds = pd.DataFrame(dict(x=x, y=y, s=size, structures=sig_str, acronym=acronym))

#cno control reversal
x = []; s = []; abbrev = []; voxels = []
for soi in sig_str:
    x.append(np.mean(df1[(df1.parent_name == soi) & (df1.Condition == "CNO_control_reversal") & (df1.Brain != "an17")]["z_score_percents"].values))
    soi = [s for s in structures if s.name==soi][0]
    if not df3.loc[df3.name == soi.name,"voxels_in_structure"].values == 0:
        voxels.append(df3.loc[df3.name == soi.name,"voxels_in_structure"].values)
    progeny = [str(xx.name) for xx in soi.progeny]
    if len(progeny) > 0:
        for progen in progeny:
            voxels.append(df3.loc[df3.name == progen,"voxels_in_structure"].values)
    s.append(np.sum(np.asarray(voxels)))
    abbrev.append(df3[df3.name == soi.name]["acronym"].values) 

size = [xx*(0.020**3)/5 for xx in s]
acronym = [iid[0] for iid in abbrev]

y = df2[df2.name.isin(sig_str)].kruskal_percent_counts_pval.values
y = np.asarray([-math.log(xx, 10) for xx in y])

df_cno_reversal = pd.DataFrame(dict(x=x, y=y, s=size, structures=sig_str, acronym=acronym))

#cno control no reversal
x = []; s = []; abbrev = []; voxels = []
for soi in sig_str:
    x.append(np.mean(df1[(df1.parent_name == soi) & (df1.Condition == "CNO_control_no_reversal")]["z_score_percents"].values))
    soi = [s for s in structures if s.name==soi][0]
    if not df3.loc[df3.name == soi.name,"voxels_in_structure"].values == 0:
        voxels.append(df3.loc[df3.name == soi.name,"voxels_in_structure"].values)
    progeny = [str(xx.name) for xx in soi.progeny]
    if len(progeny) > 0:
        for progen in progeny:
            voxels.append(df3.loc[df3.name == progen,"voxels_in_structure"].values)
    s.append(np.sum(np.asarray(voxels)))
    abbrev.append(df3[df3.name == soi.name]["acronym"].values) 

size = [xx*(0.020**3)/5 for xx in s]
acronym = [iid[0] for iid in abbrev]

y = df2[df2.name.isin(sig_str)].kruskal_percent_counts_pval.values
y = np.asarray([-math.log(xx, 10) for xx in y])

df_cno_no_reversal = pd.DataFrame(dict(x=x, y=y, s=size, structures=sig_str, acronym=acronym))

#homecage controls
x = []; s = []; abbrev = []; voxels = []
for soi in sig_str:
    x.append(np.mean(df1[(df1.parent_name == soi) & (df1.Condition == "homecage_control")]["z_score_percents"].values))
    soi = [s for s in structures if s.name==soi][0]
    if not df3.loc[df3.name == soi.name,"voxels_in_structure"].values == 0:
        voxels.append(df3.loc[df3.name == soi.name,"voxels_in_structure"].values)
    progeny = [str(xx.name) for xx in soi.progeny]
    if len(progeny) > 0:
        for progen in progeny:
            voxels.append(df3.loc[df3.name == progen,"voxels_in_structure"].values)
    s.append(np.sum(np.asarray(voxels)))
    abbrev.append(df3[df3.name == soi.name]["acronym"].values) 

size = [xx*(0.020**3)/5 for xx in s]
acronym = [iid[0] for iid in abbrev]

y = df2[df2.name.isin(sig_str)].kruskal_percent_counts_pval.values
y = np.asarray([-math.log(xx, 10) for xx in y])

df_homecage_control = pd.DataFrame(dict(x=x, y=y, s=size, structures=sig_str, acronym=acronym))
#%%
import matplotlib.patches as mpatches

fig, ax = plt.subplots(facecolor="w", figsize=(11.7, 8.3))

for key, row in df_dreadds.iterrows():
    ax.scatter(row["x"], row["y"], s=row["s"]*5, alpha=.5, color = "green", label = "DREADDs")
    ax.annotate(row["acronym"], xy=(row["x"], row["y"]), fontsize = 8, color = "midnightblue")

for key, row in df_cno_reversal.iterrows():
    ax.scatter(row["x"], row["y"], s=row["s"]*5, alpha=.5, color = "gold", label = "CNO reversal")
    ax.annotate(row["acronym"], xy=(row["x"], row["y"]), fontsize = 8, color = "darkred")

for key, row in df_cno_no_reversal.iterrows():
    ax.scatter(row["x"], row["y"], s=row["s"]*5, alpha=.5, color = "blue", label = "CNO no reversal")
    ax.annotate(row["acronym"], xy=(row["x"], row["y"]), fontsize = 8, color = "darkslategray")
    
#for key, row in df_homecage_control.iterrows():
#    ax.scatter(row["x"], row["y"], s=row["s"]*5, alpha=.5, color = "red")
#    ax.annotate(row["acronym"], xy=(row["x"], row["y"]), fontsize = 5)
#    
# Add titles (main and on axis)
plt.xlabel("Average Z-score")
plt.ylabel("-log(p-value)[Kruskal–Wallis]")
plt.title("Z-scores of Significant Structures by Kruskal–Wallis test")

green_patch = mpatches.Patch(color="green", label="DREADDs")
gold_patch = mpatches.Patch(color="gold", label="CNO reversal")
blue_patch = mpatches.Patch(color="blue", label="CNO no reversal")

plt.legend(handles=[green_patch, gold_patch, blue_patch], bbox_to_anchor=(.95, 1), loc=2, borderaxespad=0.)

plt.savefig("/home/wanglab/Desktop/bubble.svg", dpi = 300, papertype = "a3")
#%%

#PCA
from sklearn.preprocessing import StandardScaler

#make dataframe of percent counts per condition for signfiicant structures
dreadds = pd.DataFrame()
cnorev = pd.DataFrame()
cnonorev = pd.DataFrame()
homecage = pd.DataFrame()

structures = df1.parent_name.unique()

for soi in structures:
    #DREADDs
    dreadds[soi] = pd.Series(df1[(df1.parent_name == soi)& (df1.Condition == "DREADDs")].percent.values)
    dreadds["condition"] = pd.Series(df1[(df1.parent_name == soi) & (df1.Condition == "DREADDs")].Condition.values)
    #CNO reversal
    cnorev[soi] = pd.Series(df1[(df1.parent_name == soi) & (df1.Condition == "CNO_control_reversal") & (df1.Brain != "an17")].percent.values)
    cnorev["condition"] = pd.Series(df1[(df1.parent_name == soi) & (df1.Condition == "CNO_control_reversal") & (df1.Brain != "an17")].Condition.values)
    #CNO no reversal
    cnonorev[soi] = pd.Series(df1[(df1.parent_name == soi) & (df1.Condition == "CNO_control_no_reversal")].percent.values)
    cnonorev["condition"] = pd.Series(df1[(df1.parent_name == soi) & (df1.Condition == "CNO_control_no_reversal")].Condition.values)
    #homecage control
    homecage[soi] = pd.Series(df1[(df1.parent_name == soi) & (df1.Condition == "homecage_control")].percent.values)
    homecage["condition"] = pd.Series(df1[(df1.parent_name == soi) & (df1.Condition == "homecage_control")].Condition.values)
    
df_PCA = pd.concat([dreadds, cnorev, cnonorev, homecage])

#%%
# Separating out the features
arr = df_PCA.loc[:, structures].values
x = arr

# Separating out the target
y = df_PCA.loc[:,["condition"]].values

# Standardizing the features
x = StandardScaler().fit_transform(x)

#clean data
from sklearn.preprocessing import Imputer
imp = Imputer(missing_values="NaN", strategy="mean", axis=1)
cleaned_data = imp.fit_transform(x)

from sklearn.decomposition import PCA

pca = PCA(n_components=10)

principalComponents = pca.fit_transform(cleaned_data)

total_variance = np.sum(pca.explained_variance_ratio_)

finalDf = pd.concat([pd.DataFrame(principalComponents), pd.DataFrame(df_PCA["condition"].values)], axis = 1, ignore_index = True)
finalDf.to_csv("/home/wanglab/Desktop/pca_with_percents.csv")

#%%
#visualise
fig = plt.figure(figsize = (11.7,8.3))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel("Principal Component 1", fontsize = 12)
ax.set_ylabel("Principal Component 2", fontsize = 12)
ax.set_title("20 component PCA", fontsize = 12)

targets = ["DREADDs", "CNO_control_reversal", "CNO_control_no_reversal", "homecage_control"]
colors = ["g", "gold", "b", "r"]
for target, color in zip(targets,colors):
    indicesToKeep = finalDf[10] == target
    ax.scatter(finalDf.loc[indicesToKeep, 0]
               , finalDf.loc[indicesToKeep, 1]
               , c = color
               , s = 50)
ax.legend(targets)
#ax.grid()
plt.savefig("/home/wanglab/Desktop/pca.pdf")
