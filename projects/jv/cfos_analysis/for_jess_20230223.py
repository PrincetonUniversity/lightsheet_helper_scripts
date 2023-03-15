# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:08:51 2023

@author: Han
"""

import pandas as pd, pickle as pckl, numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import matplotlib as mpl
import json 

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["xtick.major.size"] = 6
mpl.rcParams["ytick.major.size"] = 6
#%%
######################################################################################################
# lobule VI
pth = r'C:/Users/Han/Desktop/lobuleVI_cellcounts.csv'
ontology_file = r'C:/Users/Han/Desktop/allen_id_table.json'
df = pd.read_csv(pth, index_col = None)

model_pth='C:/Users/Han/Desktop/model_data_contra_pma.p'

dct=pckl.load(open(model_pth  , "rb"), encoding = "latin1")
weights=dct['mat'][:,1] # only lob vi/vii (ak_pool = 1)
regions = dct['regions']

df=df[(df['Condition']=='control')]

def get_progeny(dic,parent_structure,progeny_list):
   
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

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#first calculate counts across entire nc region
av_pcounts_cfos=[];
for soi in sois:
    progeny = []
    get_progeny(ontology_dict, soi, progeny)
    dfs=[]
    for progen in progeny:
        dfs.append(df[df['name'].str.contains(progen)].groupby("Brain").sum())
    fp=pd.concat(dfs)
    fp=fp.groupby("Brain").sum() # sum up all progeny
    av_pcounts_cfos.append(np.array(fp["counts"]/df.groupby("Brain").sum()["counts"]).mean())        
av_pcounts_cfos=np.array(av_pcounts_cfos)
#sum regions per brain to get density

av_pcounts_cfos_pool = np.asarray([av_pcounts_cfos[0]+av_pcounts_cfos[1]+av_pcounts_cfos[2]+av_pcounts_cfos[3], 
        av_pcounts_cfos[5], av_pcounts_cfos[4]+av_pcounts_cfos[6], 
        av_pcounts_cfos[7]+av_pcounts_cfos[8], av_pcounts_cfos[9], av_pcounts_cfos[11], av_pcounts_cfos[10], 
        av_pcounts_cfos[12]+av_pcounts_cfos[13], av_pcounts_cfos[14]+av_pcounts_cfos[15]])


regions_cfos = ['Infralimbic, Prelimbic,\n Ant. Cingulate, Orbital',
    'Agranular insula', 'Gustatory, Visceral',
       'Somatomotor, \n Somatosensory', 'Retrosplenial', 'Visual','Post. parietal',
         'Temporal, Auditory', 'Peririhinal, Ectorhinal']
weights_wo_pp = np.delete(weights, (1)) #dropped frontal pole


#%%
# correlation between HSV model weights and DREADDs cfos 
fig, ax = plt.subplots()
ax.scatter(weights_wo_pp, np.array(av_pcounts_cfos_pool),s=80, alpha=0.3)
for i, txt in enumerate(regions_cfos):
    ax.annotate(txt, (weights_wo_pp[i], av_pcounts_cfos_pool[i]),size=8,
                xytext=(1.5,1.5),textcoords='offset points')

r,pvalr = spearmanr(weights_wo_pp, av_pcounts_cfos_pool)
p,pvalp = pearsonr(weights_wo_pp, av_pcounts_cfos_pool)
ax.annotate(f"Spearman's r={r:01f}", xy=(7,0.17), size=15)
ax.annotate(f"Spearman's p-value={pvalr:01f}", xy=(7,0.15), size=10)
ax.annotate(f"Pearson's r={p:01f}", xy=(7,0.12), size=15)
ax.annotate(f"Pearson's p-value={pvalp:01f}", xy=(7,0.10), size=10)
#find line of best fit
a, b = np.polyfit(weights_wo_pp, av_pcounts_cfos_pool, 1)
ax.plot(weights_wo_pp, a*weights_wo_pp+b, 'gray', linewidth=0.8)
ax.set_ylim([0, 0.3])
# ax.set_ylim([20,60])
plt.xlabel("Model weights for lobule VI/VII HSV injections")
plt.ylabel("Average C-fos % counts for all lobule VI")

#%%
# lobule VI thalamus
model_pth='C:/Users/Han/Desktop/thal_model_data_contra_allen.p'

dct=pckl.load(open(model_pth  , "rb"), encoding = "latin1")
weights=dct['mat'][:,1] # only lob vi/vii (ak_pool = 1)
regions = dct['regions']

#pool regions in jess dataset
sois = regions

#first calculate counts across entire nc region
av_pcounts_cfos=[];
for soi in sois:
    progeny = []
    get_progeny(ontology_dict, soi, progeny)
    dfs=[]
    progeny.append(soi) # to account for thalamic nuclei with no children        
    for progen in progeny:
        dfs.append(df[df['name'].str.contains(progen)].groupby("Brain").sum())
    try:
        fp=pd.concat(dfs)
        fp=fp.groupby("Brain").sum() # sum up all progeny
        av_pcounts_cfos.append(np.nanmean(np.array(fp["counts"]/df.groupby("Brain").sum()["counts"])))
        # nanmean to account for brains with no thalamic nuclei counts
    except:
        print(progen)
av_pcounts_cfos=np.array(av_pcounts_cfos)

#abbreviations
regions_cfos = np.array(["VPM", "VPL","VA-L", "AV",  "LD",
       "PV", "MedHab", "LP",
       "PoT", "MD", "Po",
       "VM", "RTN"])

#%%
# correlation between HSV model weights and DREADDs cfos 
fig, ax = plt.subplots()
ax.scatter(weights, np.array(av_pcounts_cfos),s=80, alpha=0.3)
for i, txt in enumerate(regions_cfos):
    ax.annotate(txt, (weights[i], av_pcounts_cfos[i]),size=8,
                xytext=(-1.5,1.5),textcoords='offset points',
                horizontalalignment='right',
            verticalalignment='bottom')

r,pvalr = spearmanr(weights, av_pcounts_cfos)
p,pvalp = pearsonr(weights, av_pcounts_cfos)
ax.annotate(f"Spearman's r={r:01f}", xy=(4.5,0.008), size=15)
ax.annotate(f"Spearman's p-value={pvalr:01f}", xy=(4.5,0.007), size=10)
ax.annotate(f"Pearson's r={p:01f}", xy=(4.5,0.006), size=15)
ax.annotate(f"Pearson's p-value={pvalp:01f}", xy=(4.5,0.005), size=10)
#find line of best fit
a, b = np.polyfit(weights, av_pcounts_cfos, 1)
ax.plot(weights, a*weights+b, 'gray', linewidth=0.8)

# ax.set_ylim([20,60])
plt.xlabel("Model weights for lobule VI/VII HSV injections")
plt.ylabel("Average C-fos % counts for lobule VI controls")

# %%
######################################################################################################
# crus I neocortex
pth = r'C:/Users/Han/Desktop/cell_counts_dataframe_w_percents_density_somecrusI.csv'
df = pd.read_csv(pth, index_col = None)

model_pth='C:/Users/Han/Desktop/model_data_contra_pma.p'

dct=pckl.load(open(model_pth  , "rb"), encoding = "latin1")
weights=dct['mat'][:,1] # only lob vi/vii (ak_pool = 1)
regions = dct['regions']

df=df[(df['Condition']=='DREADDs')]

def get_progeny(dic,parent_structure,progeny_list):
   
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

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#first calculate counts across entire nc region
av_pcounts_cfos=[];
for soi in sois:
    progeny = []
    get_progeny(ontology_dict, soi, progeny)
    dfs=[]
    for progen in progeny:
        dfs.append(df[df['name'].str.contains(progen)].groupby("Brain").sum())
    fp=pd.concat(dfs)
    fp=fp.groupby("Brain").sum() # sum up all progeny
    av_pcounts_cfos.append(np.nanmean(np.array(fp["counts"]/df.groupby("Brain").sum()["counts"])))        
av_pcounts_cfos=np.array(av_pcounts_cfos)
#sum regions per brain to get density

av_pcounts_cfos_pool = np.asarray([av_pcounts_cfos[0]+av_pcounts_cfos[1]+av_pcounts_cfos[2]+av_pcounts_cfos[3], 
        av_pcounts_cfos[5], av_pcounts_cfos[4]+av_pcounts_cfos[6], 
        av_pcounts_cfos[7]+av_pcounts_cfos[8], av_pcounts_cfos[9], av_pcounts_cfos[11], av_pcounts_cfos[10], 
        av_pcounts_cfos[12]+av_pcounts_cfos[13], av_pcounts_cfos[14]+av_pcounts_cfos[15]])


regions_cfos = ['Infralimbic, Prelimbic,\n Ant. Cingulate, Orbital',
    'Agranular insula', 'Gustatory, Visceral',
       'Somatomotor, \n Somatosensory', 'Retrosplenial', 'Visual','Post. parietal',
         'Temporal, Auditory', 'Peririhinal, Ectorhinal']
weights_wo_pp = np.delete(weights, (1)) #dropped frontal pole


#%%
fig, ax = plt.subplots()
ax.scatter(weights_wo_pp, np.array(av_pcounts_cfos_pool),s=80, alpha=0.3)
for i, txt in enumerate(regions_cfos):
    ax.annotate(txt, (weights_wo_pp[i], av_pcounts_cfos_pool[i]),size=8,
                xytext=(1.5,1.5),textcoords='offset points')

r,pvalr = spearmanr(weights_wo_pp, av_pcounts_cfos_pool)
p,pvalp = pearsonr(weights_wo_pp, av_pcounts_cfos_pool)
ax.annotate(f"Spearman's r={r:01f}", xy=(7,0.17), size=15)
ax.annotate(f"Spearman's p-value={pvalr:01f}", xy=(7,0.15), size=10)
ax.annotate(f"Pearson's r={p:01f}", xy=(7,0.12), size=15)
ax.annotate(f"Pearson's p-value={pvalp:01f}", xy=(7,0.10), size=10)
#find line of best fit
a, b = np.polyfit(weights_wo_pp, av_pcounts_cfos_pool, 1)
ax.plot(weights_wo_pp, a*weights_wo_pp+b, 'gray', linewidth=0.8)
ax.set_ylim([0, 0.3])
plt.xlabel("Model weights for Crus I/II HSV injections")
plt.ylabel("Average C-fos % counts for Crus I CN0 reversal controls")

#%%
# crus I thalamus
model_pth='C:/Users/Han/Desktop/thal_model_data_contra_allen.p'

dct=pckl.load(open(model_pth  , "rb"), encoding = "latin1")
weights=dct['mat'][:,1] # only lob vi/vii (ak_pool = 1)
regions = dct['regions']

#pool regions in jess dataset
sois = regions

#first calculate counts across entire nc region
av_pcounts_cfos=[];
for soi in sois:
    progeny = []
    get_progeny(ontology_dict, soi, progeny)
    dfs=[]
    progeny.append(soi) # to account for thalamic nuclei with no children        
    for progen in progeny:
        dfs.append(df[df['name'].str.contains(progen)].groupby("Brain").sum())
    try:
        fp=pd.concat(dfs)
        fp=fp.groupby("Brain").sum() # sum up all progeny
        av_pcounts_cfos.append(np.nanmean(np.array(fp["counts"]/df.groupby("Brain").sum()["counts"])))
        # nanmean to account for brains with no thalamic nuclei counts
    except:
        print(progen)
av_pcounts_cfos=np.array(av_pcounts_cfos)

#abbreviations
regions_cfos = np.array(["VPM", "VPL","VA-L", "AV",  "LD",
       "PV", "MedHab", "LP",
       "PoT", "MD", "Po",
       "VM", "RTN"])

#%%
# correlation between HSV model weights and DREADDs cfos 
fig, ax = plt.subplots()
ax.scatter(weights, np.array(av_pcounts_cfos),s=80, alpha=0.3)
for i, txt in enumerate(regions_cfos):
    ax.annotate(txt, (weights[i], av_pcounts_cfos[i]),size=8,
                xytext=(-1.5,1.5),textcoords='offset points',
                horizontalalignment='right',
            verticalalignment='bottom')

r,pvalr = spearmanr(weights, av_pcounts_cfos)
p,pvalp = pearsonr(weights, av_pcounts_cfos)
ax.annotate(f"Spearman's r={r:01f}", xy=(4.5,0.008), size=15)
ax.annotate(f"Spearman's p-value={pvalr:01f}", xy=(4.5,0.007), size=10)
ax.annotate(f"Pearson's r={p:01f}", xy=(4.5,0.006), size=15)
ax.annotate(f"Pearson's p-value={pvalp:01f}", xy=(4.5,0.005), size=10)
#find line of best fit
a, b = np.polyfit(weights, av_pcounts_cfos, 1)
ax.plot(weights, a*weights+b, 'gray', linewidth=0.8)

ax.set_ylim([0,0.01])
plt.xlabel("Model weights for crus I/II HSV injections")
plt.ylabel("Average C-fos % counts for Crus I DREADDs")
