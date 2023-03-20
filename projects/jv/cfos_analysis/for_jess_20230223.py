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
pth = r'C:/Users/Han/Desktop/jess_df.xlsx'
ontology_file = r'C:/Users/Han/Desktop/allen_id_table.json'
df = pd.read_excel(pth, index_col = None)

model_pth='C:/Users/Han/Desktop/model_data_contra_pma.p'

dct=pckl.load(open(model_pth  , "rb"), encoding = "latin1")
weights=dct['mat'][:,1] # only lob vi/vii (ak_pool = 1)
regions = dct['regions']

sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Orbital area", 
            "Gustatory areas", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]

#per condition
dfgroup=df.groupby("condition").mean()[sois]
dfgroup['Infralimbic, Prelimbic,\n Ant. Cingulate, Orbital']=dfgroup[["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Orbital area"]].sum(axis=1)
dfgroup['Gustatory, Visceral']=dfgroup[["Gustatory areas", "Visceral area"]].sum(axis=1)
dfgroup['Somatomotor, Somatosensory']=dfgroup[["Somatosensory areas", "Somatomotor areas"]].sum(axis=1)
dfgroup['Tenporal, Auditory']=dfgroup[["Temporal association areas","Auditory areas"]].sum(axis=1)
dfgroup['Peririhinal, Ectorhinal']=dfgroup[["Ectorhinal area", "Perirhinal area"]].sum(axis=1)
#sum regions per brain to get density

regions_cfos = ['Infralimbic, Prelimbic,\n Ant. Cingulate, Orbital','Gustatory, Visceral',
       'Somatomotor, Somatosensory', 'Retrosplenial', 'Visual','Post. parietal',
         'Temporal, Auditory', 'Peririhinal, Ectorhinal']
weights_wo_pp = np.delete(weights, (1,2)) #dropped frontal pole and insula
av_pcounts = dfgroup[['Infralimbic, Prelimbic,\n Ant. Cingulate, Orbital','Gustatory, Visceral',
       'Somatomotor, Somatosensory',"Retrosplenial area", 
        "Visual areas", "Posterior parietal association areas",  'Tenporal, Auditory',
        'Peririhinal, Ectorhinal']]

#%%
# correlation between HSV model weights and DREADDs cfos 
weights=dct['mat'][:,1] # only lob vi/vii (ak_pool = 1) or crus I (ak_pool=-2)
weights=dct['mat'][:,-2] # only lob vi/vii (ak_pool = 1) or crus I (ak_pool=-2)
weights_wo_pp = np.delete(weights, (1,2)) #dropped frontal pole and insula

fig, ax = plt.subplots()
Y = av_pcounts[av_pcounts.index=='CNO_control_reversal'].values[0] #change condition

ax.scatter(weights_wo_pp, Y,
           s=80, alpha=0.3)
for i, txt in enumerate(regions_cfos):
    ax.annotate(txt, (weights_wo_pp[i], Y[i]),size=8,
                xytext=(1.5,1.5),textcoords='offset points')

r,pvalr = spearmanr(weights_wo_pp, Y)
p,pvalp = pearsonr(weights_wo_pp, Y)
ax.annotate(f"Spearman's r={r:01f}", xy=(7,0.17), size=15)
ax.annotate(f"Spearman's p-value={pvalr:01f}", xy=(7,0.15), size=10)
ax.annotate(f"Pearson's r={p:01f}", xy=(7,0.12), size=15)
ax.annotate(f"Pearson's p-value={pvalp:01f}", xy=(7,0.10), size=10)
#find line of best fit
a, b = np.polyfit(weights_wo_pp, Y, 1)
ax.plot(weights_wo_pp, a*weights_wo_pp+b, 'gray', linewidth=0.8)
# ax.set_ylim([20,60])
plt.xlabel("Model weights for crus I/II HSV injections") #lobule VI/VII
plt.ylabel("Average C-fos % counts for crus I") #lobule VI
# plt.xlabel("Model weights for lobule VI/VII HSV injections") #
# plt.ylabel("Average C-fos % counts for lobule VI") #


#%%
# thalamus
model_pth='C:/Users/Han/Desktop/thal_model_data_contra_allen.p'

dct=pckl.load(open(model_pth  , "rb"), encoding = "latin1")
regions = dct['regions']

#pool regions in jess dataset
sois = ['Ventral posteromedial nucleus of the thalamus',
       'Ventral posterolateral nucleus of the thalamus',
       'Ventral anterior-lateral complex of the thalamus',
       'Anteroventral nucleus of thalamus',
       'Lateral dorsal nucleus of thalamus',
       'Paraventricular nucleus of the thalamus', 'Medial habenula',
       'Lateral posterior nucleus of the thalamus',
       'Posterior triangular thalamic nucleus',
       'Posterior complex of the thalamus',
       'Ventral medial nucleus of the thalamus',
       'Reticular nucleus of the thalamus']
av_pcounts=df.groupby("condition").mean()[sois]
#abbreviations
regions_cfos = np.array(["VPM", "VPL","VA-L", "AV",  "LD",
       "PV", "MedHab", "LP",
       "PoT", "Po",
       "VM", "RTN"])

#%%
# plot
weights=dct['mat'][:,1] # only lob vi/vii (ak_pool = 1)
weights_wo_pp = np.delete(weights, 4) # dropped md

fig, ax = plt.subplots()
Y = av_pcounts[av_pcounts.index=='CNO_control_reversal'].values[0]

ax.scatter(weights_wo_pp, Y,
           s=80, alpha=0.3)
for i, txt in enumerate(regions_cfos):
    ax.annotate(txt, (weights_wo_pp[i], Y[i]),size=8,
                xytext=(1.5,1.5),textcoords='offset points')

r,pvalr = spearmanr(weights_wo_pp, Y)
p,pvalp = pearsonr(weights_wo_pp, Y)
ax.annotate(f"Spearman's r={r:01f}", xy=(4,0.003), size=15)
ax.annotate(f"Spearman's p-value={pvalr:01f}", xy=(4,0.0025), size=10)
ax.annotate(f"Pearson's r={p:01f}", xy=(4,0.002), size=15)
ax.annotate(f"Pearson's p-value={pvalp:01f}", xy=(4,0.0015), size=10)
#find line of best fit
a, b = np.polyfit(weights_wo_pp, Y, 1)
ax.plot(weights_wo_pp, a*weights_wo_pp+b, 'gray', linewidth=0.8)
# ax.set_ylim([0,0.01])
# plt.xlabel("Model weights for crus I/II HSV injections") #lobule VI/VII
# plt.ylabel("Average C-fos % counts for crus I") #lobule VI
plt.xlabel("Model weights for lobule VI/VII HSV injections") #
plt.ylabel("Average C-fos % counts for lobule VI") #

