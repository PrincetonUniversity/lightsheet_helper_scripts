# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 12:50:55 2020

@author: Zahra
"""

import os, pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, json

dst = r"Z:\zahra\h129_contra_vs_ipsi"
df_pth = r"Y:\atlas\ls_id_table_w_voxelcounts.xlsx"
cells_regions_pth = os.path.join(dst, "data/thal_contra_counts_23_brains_80um_ventric_erosion.csv")

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.025
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

brains = cells_regions.columns[:-1]

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
ontology_file = r"Y:\atlas\allen_atlas\allen.json"

with open(ontology_file) as json_file:
    ontology_dict = json.load(json_file)

sois = ["Thalamus", 
         "Ventral anterior-lateral complex of the thalamus",
         "Ventral medial nucleus of the thalamus",
         "Ventral posterolateral nucleus of the thalamus",
         "Ventral posteromedial nucleus of the thalamus",
         "Posterior triangular thalamic nucleus",
         "Subparafascicular nucleus",
         "Subparafascicular area",
         "Peripeduncular nucleus",
         "Medial geniculate complex",
         "Dorsal part of the lateral geniculate complex",
         "Lateral posterior nucleus of the thalamus",
         "Posterior complex of the thalamus",
         "Posterior limiting nucleus of the thalamus",
         "Suprageniculate nucleus",
         "Ethmoid nucleus of the thalamus",
         "Retroethmoid nucleus",
         "Anteroventral nucleus of thalamus",
         "Anteromedial nucleus",
         "Anterodorsal nucleus",
         "Interanteromedial nucleus of the thalamus",
         "Interanterodorsal nucleus of the thalamus",
         "Lateral dorsal nucleus of thalamus",
         "Intermediodorsal nucleus of the thalamus",
         "Mediodorsal nucleus of thalamus",
         "Submedial nucleus of the thalamus",
         "Perireunensis nucleus",
         "Paraventricular nucleus of the thalamus",
         "Parataenial nucleus",
         "Nucleus of reuniens",
         "Xiphoid thalamic nucleus",
         "Rhomboid nucleus",
         "Central medial nucleus of the thalamus",
         "Paracentral nucleus",
         "Central lateral nucleus of the thalamus",
         "Parafascicular nucleus",
         "Posterior intralaminar thalamic nucleus",
         "Reticular nucleus of the thalamus",
         "Intergeniculate leaflet of the lateral geniculate complex",
         "Intermediate geniculate nucleus",
         "Ventral part of the lateral geniculate complex",
         "Subgeniculate nucleus",
         "Medial habenula",
         "Lateral habenula",
         "Pineal body"]

#first calculate counts across entire region
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

#voxels
vol = []
for soi in sois[1:]:
    progeny = []; counts = []; iids = []
    get_progeny(ontology_dict, soi, progeny)
    #add counts from main structure
    try:
        counts.append(ann_df.loc[ann_df.name == soi, "voxels_in_structure"].values[0])
    except:
        print(soi)
    for progen in progeny:
        counts.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    vol.append(np.array(counts).sum(axis = 0))
vol = np.array(vol)        

density = np.nan_to_num(np.array([xx/(vol[i]*(scale_factor**3)) for i, xx in enumerate(counts_per_struct[1:])]).T) #remove thalamus

pcounts = np.nan_to_num(np.asarray([((brain[1:]/brain[0])*100) for brain in counts_per_struct.T]))    

#%%
#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
order = np.argsort(np.median(pcounts, axis = 0))[::-1]
sois_sort = np.array(sois[1:])[order][:20]

#boxplots of percent counts
plt.figure(figsize = (5,7))
df = pd.DataFrame(pcounts)
df.columns = sois[1:]
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort)
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={'facecolor':'None'}, order = sois_sort)
plt.xlabel("% of total thalamic cells")
plt.ylabel("Thalamic nuclei")
plt.savefig(os.path.join(dst, "thal_pcounts_boxplots.pdf"), bbox_inches = "tight")

order = np.argsort(np.median(density, axis = 0))[::-1]
sois_sort = np.array(sois[1:])[order][:20]

#boxplots of percent counts
plt.figure(figsize = (5,7))
df = pd.DataFrame(density)
df.columns = sois[1:]
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort)
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={'facecolor':'None'}, order = sois_sort)
plt.xlabel("Cells/$mm^3$")
plt.ylabel("Thalamic nuclei")
plt.savefig(os.path.join(dst, "thal_density_boxplots.pdf"), bbox_inches = "tight")

#%%
#save out for hjb

df = pd.DataFrame(pcounts.T)
df.columns = brains
df["Structure"] = sois[1:]

df.to_csv(os.path.join(dst, "pcounts_thalamus_contralateral_hsv_80um_erosion.csv"))

df = pd.DataFrame(density.T)
df.columns = brains
df["Structure"] = sois[1:]
df["Voxels in structure"] = vol

df.to_csv(os.path.join(dst, "density_thalamus_contralateral_hsv_80um_erosion.csv"))