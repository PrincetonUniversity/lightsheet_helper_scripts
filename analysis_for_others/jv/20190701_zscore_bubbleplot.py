#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 19:25:35 2019

@author: wanglab
"""

#Bubbleplot
import math, pandas as pd, numpy as np, matplotlib as mpl, os
import matplotlib.pyplot as plt
from lightsheet.network_analysis import make_structure_objects
from scipy.stats import zscore

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#set paths here
lst_pth = "/jukebox/LightSheetTransfer/atlas/"
data_pth = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/more_selected_structures/"
dst = "/home/wanglab/Desktop"

#no need to run more than once
#get structure ontology
structures = make_structure_objects(os.path.join(lst_pth, "ls_id_table_w_voxelcounts.xlsx"), 
                                    remove_childless_structures_not_repsented_in_ABA = True, 
                                    ann_pth=os.path.join(lst_pth, "annotation_sagittal_atlas_20um_iso.tif"))

#%%
#bubble plot
df1 = pd.read_csv(os.path.join(data_pth, "select_structures_percent_counts_for_plots.csv"))
df2 = pd.read_csv(os.path.join(data_pth, "one_way_anova_all_structures.csv"))
df3 = pd.read_excel(os.path.join(lst_pth, "ls_id_table_w_voxelcounts.xlsx"))
ann_pth = os.path.join(lst_pth, "annotation_sagittal_atlas_20um_iso.tif")

#drop homecage ctrl
df1 = df1[df1.condition != "homecage_control"]
#look at stuff per condition
sig_str = df2[(df2.anova_percent_counts_pval < 0.1)].name.values

#calc z score
structs = df1.name.unique()

for nm in structs:
    df1.loc[(df1.name == nm), "z_score_percents"] = zscore(df1[df1.name == nm]["percent"])

#how big bubbles are
sizef = 3

#dreadds
x = []; s = []; abbrev = []; voxels = []
for soi in sig_str:
    x.append(np.mean(df1[(df1.name == soi) & (df1.condition == "DREADDs")]["z_score_percents"].values))
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

size = [xx*(0.020**3)/sizef for xx in s]
acronym = [iid[0] for iid in abbrev]

y = df2[df2.name.isin(sig_str)].anova_percent_counts_pval.values
y = np.asarray([-math.log(xx, 10) for xx in y])

df_dreadds = pd.DataFrame(dict(x=x, y=y, s=size, structures=sig_str, acronym=acronym))

#cno control reversal
x = []; s = []; abbrev = []; voxels = []
for soi in sig_str:
    x.append(np.mean(df1[(df1.name == soi) & (df1.condition == "CNO_control_reversal") & (df1.animal != "an17")]["z_score_percents"].values))
    soi = [s for s in structures if s.name==soi][0]
    if not df3.loc[df3.name == soi.name,"voxels_in_structure"].values == 0:
        voxels.append(df3.loc[df3.name == soi.name,"voxels_in_structure"].values)
    progeny = [str(xx.name) for xx in soi.progeny]
    if len(progeny) > 0:
        for progen in progeny:
            voxels.append(df3.loc[df3.name == progen,"voxels_in_structure"].values)
    s.append(np.sum(np.asarray(voxels)))
    abbrev.append(df3[df3.name == soi.name]["acronym"].values) 

size = [xx*(0.020**3)/sizef for xx in s]
acronym = [iid[0] for iid in abbrev]

y = df2[df2.name.isin(sig_str)].anova_percent_counts_pval.values
y = np.asarray([-math.log(xx, 10) for xx in y])

df_cno_reversal = pd.DataFrame(dict(x=x, y=y, s=size, structures=sig_str, acronym=acronym))

#cno control no reversal
x = []; s = []; abbrev = []; voxels = []
for soi in sig_str:
    x.append(np.mean(df1[(df1.name == soi) & (df1.condition == "CNO_control_no_reversal")]["z_score_percents"].values))
    soi = [s for s in structures if s.name==soi][0]
    if not df3.loc[df3.name == soi.name,"voxels_in_structure"].values == 0:
        voxels.append(df3.loc[df3.name == soi.name,"voxels_in_structure"].values)
    progeny = [str(xx.name) for xx in soi.progeny]
    if len(progeny) > 0:
        for progen in progeny:
            voxels.append(df3.loc[df3.name == progen,"voxels_in_structure"].values)
    s.append(np.sum(np.asarray(voxels)))
    abbrev.append(df3[df3.name == soi.name]["acronym"].values) 

size = [xx*(0.020**3)/sizef for xx in s]
acronym = [iid[0] for iid in abbrev]

y = df2[df2.name.isin(sig_str)].anova_percent_counts_pval.values
y = np.asarray([-math.log(xx, 10) for xx in y])

df_cno_no_reversal = pd.DataFrame(dict(x=x, y=y, s=size, structures=sig_str, acronym=acronym))


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
plt.ylabel("-log(p-value)[ANOVA]")
plt.title("One-way ANOVA Significant Structures")

green_patch = mpatches.Patch(color="green", label="DREADDs")
gold_patch = mpatches.Patch(color="gold", label="CNO reversal")
blue_patch = mpatches.Patch(color="blue", label="CNO no reversal")

plt.hlines(y=-np.log10(0.05), colors = 'grey', linestyles = "dashed", xmin=-1, xmax=1)
plt.hlines(y=-np.log10(0.01), colors = 'k', linestyles = "dotted", xmin=-1, xmax=1)

plt.legend(handles=[green_patch, gold_patch, blue_patch], bbox_to_anchor=(.65, 1), loc=2, borderaxespad=0.)

plt.savefig(os.path.join(dst, "bubble.pdf"), dpi = 300, papertype = "a3")