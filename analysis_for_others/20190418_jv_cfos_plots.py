# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 22:09:13 2019

@author: Zahra
"""

import pandas as pd, seaborn as sns, matplotlib.pyplot as plt
import os
import numpy as np

sum_pth = '/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/summed_cell_counts_dataframe_w_zscores_per_structure.csv'
   
df = pd.read_csv(sum_pth)

sns.set(font_scale=0.3, style = "white")#keep font small

#shared x and y axis
fig, ax = plt.subplots(1,4, sharex = True, sharey = True)

#boxplot z score
plt.subplot(141)
ax = sns.boxplot(x="z_score_percents", y="parent", data=df[df.Condition == "homecage_control"], showfliers = False, color = "r")
ax = sns.swarmplot(x="z_score_percents", y="parent", data=df[df.Condition == "homecage_control"], color="k", size = 2)

plt.subplot(142)
#boxplot z score
ax = sns.boxplot(x="z_score_percents", y="parent", data=df[df.Condition == "CNO_control_no_reversal"], showfliers = False, color = "royalblue")
ax = sns.swarmplot(x="z_score_percents", y="parent", data=df[df.Condition == "CNO_control_no_reversal"], color="k", size = 2)

ax.get_yaxis().set_visible(False)

plt.subplot(143)
#boxplot z score
ax = sns.boxplot(x="z_score_percents", y="parent", data=df[df.Condition == "CNO_control_reversal"], showfliers = False, color = "sandybrown")
ax = sns.swarmplot(x="z_score_percents", y="parent", data=df[df.Condition == "CNO_control_reversal"], color="k", size = 2)

ax.get_yaxis().set_visible(False)

plt.subplot(144)
#boxplot z score
ax = sns.boxplot(x="z_score_percents", y="parent", data=df[df.Condition == "DREADDs"], showfliers = False, color = "green")
ax = sns.swarmplot(x="z_score_percents", y="parent", data=df[df.Condition == "DREADDs"], color="k", size = 2)

ax.get_yaxis().set_visible(False)

plt.savefig("/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/z_score_boxplots/summed_structures_percent_per_total_count_zscores.pdf", bbox_inches = "tight")   

#%%
sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area", "Entorhinal area"]

#only looking at neocortical structures now
sns.set(font_scale=0.7, style = "white")#keep font small

all_pth = '/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/z_score_boxplots/csv_files/z_score_all.csv'

df = pd.read_csv(all_pth, index_col = None)

def plot_specified_sois(df, sois, hierarchy, dst):
     
    """ 
    plots boxplots side by side 
    inputs: 
        df = DREADDs data frame
        sois = allen structures
        hierarchy = neocortex, thalamus, etc (STRING)
        dst = destination directory
    """
    
    #shared x and y axis
    fig, ax = plt.subplots(1,4, sharex = True, sharey = True)
    
    #boxplot z score
    plt.subplot(141)
    ax = sns.boxplot(x="z_score_percents", y="name", data=df[(df.Condition == "homecage_control") & (df.name.isin(sois))], showfliers = False, color = "r")
    ax = sns.swarmplot(x="z_score_percents", y="name", data=df[(df.Condition == "homecage_control") & (df.name.isin(sois))], color="k", size = 2)
    
    #formatting
    plt.axvline(0, color = 'grey')
    
    plt.subplot(142)
    ax = sns.boxplot(x="z_score_percents", y="name", data=df[(df.Condition == "CNO_control_no_reversal") & (df.name.isin(sois))], showfliers = False, color = "royalblue")
    ax = sns.swarmplot(x="z_score_percents", y="name", data=df[(df.Condition == "CNO_control_no_reversal") & (df.name.isin(sois))], color="k", size = 2)
    
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('')
    plt.axvline(0, color = 'grey')
    
    plt.subplot(143)
    ax = sns.boxplot(x="z_score_percents", y="name", data=df[(df.Condition == "CNO_control_reversal") & (df.name.isin(sois))], showfliers = False, color = "sandybrown")
    ax = sns.swarmplot(x="z_score_percents", y="name", data=df[(df.Condition == "CNO_control_reversal") & (df.name.isin(sois))], color="k", size = 2)
    
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('')
    plt.axvline(0, color = 'grey')
    
    plt.subplot(144)
    ax = sns.boxplot(x="z_score_percents", y="name", data=df[(df.Condition == "DREADDs") & (df.name.isin(sois))], showfliers = False, color = "green")
    ax = sns.swarmplot(x="z_score_percents", y="name", data=df[(df.Condition == "DREADDs") & (df.name.isin(sois))], color="k", size = 2)
    
    #formatting
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('')
    plt.axvline(0, color = 'grey')
    
    plt.savefig(os.path.join(dst, "summed_structures_boxplot_percent_per_total_count_zscores_{}.pdf".format(hierarchy)), bbox_inches = "tight")   
    
    #look at density
    #shared x and y axis
    fig, ax = plt.subplots(1,4, sharex = True, sharey = True)
    
    #boxplot z score
    plt.subplot(141)
    ax = sns.boxplot(x="z_score_density", y="name", data=df[(df.Condition == "homecage_control") & (df.name.isin(sois))], showfliers = False, color = "r")
    ax = sns.swarmplot(x="z_score_density", y="name", data=df[(df.Condition == "homecage_control") & (df.name.isin(sois))], color="k", size = 2)

    #formatting
    plt.axvline(0, color = 'grey')
    
    plt.subplot(142)
    ax = sns.boxplot(x="z_score_density", y="name", data=df[(df.Condition == "CNO_control_no_reversal") & (df.name.isin(sois))], showfliers = False, color = "royalblue")
    ax = sns.swarmplot(x="z_score_density", y="name", data=df[(df.Condition == "CNO_control_no_reversal") & (df.name.isin(sois))], color="k", size = 2)
    
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('')
    plt.axvline(0, color = 'grey')
    
    plt.subplot(143)
    ax = sns.boxplot(x="z_score_density", y="name", data=df[(df.Condition == "CNO_control_reversal") & (df.name.isin(sois))], showfliers = False, color = "sandybrown")
    ax = sns.swarmplot(x="z_score_density", y="name", data=df[(df.Condition == "CNO_control_reversal") & (df.name.isin(sois))], color="k", size = 2)
    
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('')
    plt.axvline(0, color = 'grey')
    
    plt.subplot(144)
    ax = sns.boxplot(x="z_score_density", y="name", data=df[(df.Condition == "DREADDs") & (df.name.isin(sois))], showfliers = False, color = "green")
    ax = sns.swarmplot(x="z_score_density", y="name", data=df[(df.Condition == "DREADDs") & (df.name.isin(sois))], color="k", size = 2)
    
    #formatting
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('')
    plt.axvline(0, color = 'grey')
    
    plt.savefig(os.path.join(dst, "summed_structures_boxplot_density_zscores_{}.pdf".format(hierarchy)), bbox_inches = "tight")

plot_specified_sois(df, sois, "neocortex", "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/z_score_boxplots")

#now looking at thalamus, etc.
sois = ["Ventral posterior complex of the thalamus", "Ventral posterolateral nucleus of the thalamus",
        "Ventral posteromedial nucleus of the thalamus",
            "Ventral anterior-lateral complex of the thalamus", "Medial geniculate complex",
            "Lateral group of the dorsal thalamus", "Anterior group of the dorsal thalamus", "Medial group of the dorsal thalamus",
            "Intralaminar nuclei of the dorsal thalamus", "Geniculate group, ventral thalamus", "Ventral group of the dorsal thalamus",
            "Reticular nucleus of the thalamus"]

plot_specified_sois(df, sois, "thalamus", "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/z_score_boxplots")

#amygdala
sois = ["Basolateral amygdalar nucleus", "Basomedial amygdalar nucleus", "Posterior amygdalar nucleus", "Lateral amygdalar nucleus"]

plot_specified_sois(df, sois, "amgdala", "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/z_score_boxplots")

#hypothalamic nuclei
sois = ["Anterior hypothalamic nucleus", "Ventromedial hypothalamic nucleus", "Posterior hypothalamic nucleus",
        "Lateral hypothalamic area"]

plot_specified_sois(df, sois, "hypothalamus", "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/z_score_boxplots")

#%%

import math

#bubble plot
df1 = pd.read_csv('/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/summed_cell_counts_dataframe_w_zscores_per_structure.csv')
df2 = pd.read_csv('/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/two_way_anova_pooled_structures_reversal_vs_no_reversal.csv')

#look at DREADDs
sig_str = df2[df2.two_way_pval < 0.1].parent.values

x = []; s = []; ids = []
for sig in sig_str:
    x.append(np.mean(df1[(df1.parent == sig) & (df1.Condition == "DREADDs")]["z_score_percents"].values))
    s.append(df1[(df1.parent == sig) & (df1.Condition == "DREADDs")]["volume"].drop_duplicates().values)
    ids.append(df1[(df1.parent == sig) & (df1.Condition == "DREADDs")]["id"].drop_duplicates().values) 

#formatting
s = [xx*30 for xx in s]
ids = [iid[0] for iid in ids]

y = df2[df2.parent.isin(sig_str)].two_way_pval.values
y = np.asarray([-math.log(xx, 10) for xx in y])

df = pd.DataFrame(dict(x=x, y=y, s=s, structures=sig_str, ids=ids))

fig, ax = plt.subplots(facecolor='w')

for key, row in df.iterrows():
    ax.scatter(row['x'], row['y'], s=row['s']*5, alpha=.5)
    ax.annotate(row['structures'], xy=(row['x'], row['y']), fontsize = 5)
 
# Add titles (main and on axis)
plt.xlabel("Average Z-score")
plt.ylabel("-log(p-value)[ANOVA]")
plt.title("DREADDs")
 
plt.savefig("/home/wanglab/Desktop/test.pdf", dpi = 300)
