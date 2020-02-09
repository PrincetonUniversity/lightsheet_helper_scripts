# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 12:50:55 2020

@author: Zahra
"""

import os, pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, json, matplotlib as mpl

#TP
plt.rcParams["axes.grid"] = False

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
cells_regions_pth = os.path.join(dst, "data/thal_contra_counts_23_brains_80um_ventric_erosion.csv")

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
scale_factor = 0.025
ann_df = pd.read_excel(df_pth).drop(columns = ["Unnamed: 0"])

brains = cells_regions.columns[:-1]

def get_progeny(dic,parent_structure,progeny_list):
   
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

sois_dict = {"Thalamus": "Sensory-motor",
         "Ventral anterior-lateral complex of the thalamus": "Sensory-motor", 
         "Ventral medial nucleus of the thalamus": "Sensory-motor", 
         "Ventral posterolateral nucleus of the thalamus": "Sensory-motor", 
         "Ventral posteromedial nucleus of the thalamus": "Sensory-motor", 
         "Posterior triangular thalamic nucleus": "Sensory-motor", 
         "Peripeduncular nucleus": "Sensory-motor", 
         "Medial geniculate complex": "Sensory-motor",
         "Dorsal part of the lateral geniculate complex": "Sensory-motor", 
         "Lateral posterior nucleus of the thalamus": "Polymodal association", 
         "Posterior complex of the thalamus": "Polymodal association", 
         "Suprageniculate nucleus": "Polymodal association", 
         "Ethmoid nucleus of the thalamus": "Polymodal association", 
         "Retroethmoid nucleus": "Polymodal association", 
         "Anteroventral nucleus of thalamus": "Polymodal association", 
         "Anteromedial nucleus": "Polymodal association", 
         "Anterodorsal nucleus": "Polymodal association", 
         "Interanteromedial nucleus of the thalamus": "Polymodal association", 
         "Interanterodorsal nucleus of the thalamus": "Polymodal association", 
         "Lateral dorsal nucleus of thalamus": "Polymodal association", 
         "Intermediodorsal nucleus of the thalamus": "Polymodal association", 
         "Mediodorsal nucleus of thalamus": "Polymodal association", 
         "Submedial nucleus of the thalamus": "Polymodal association", 
         "Perireunensis nucleus": "Polymodal association", 
         "Paraventricular nucleus of the thalamus": "Polymodal association", 
         "Parataenial nucleus": "Polymodal association", 
         "Nucleus of reuniens": "Polymodal association", 
         "Xiphoid thalamic nucleus": "Polymodal association", 
         "Rhomboid nucleus": "Polymodal association", 
         "Central medial nucleus of the thalamus": "Polymodal association", 
         "Paracentral nucleus": "Polymodal association", 
         "Central lateral nucleus of the thalamus": "Polymodal association", 
         "Parafascicular nucleus": "Polymodal association", 
         "Posterior intralaminar thalamic nucleus": "Polymodal association", 
         "Reticular nucleus of the thalamus": "Polymodal association", 
         "Intermediate geniculate nucleus": "Polymodal association", 
         "Ventral part of the lateral geniculate complex": "Polymodal association", 
         "Subgeniculate nucleus": "Polymodal association", 
         "Medial habenula": "Polymodal association", 
         "Lateral habenula": "Polymodal association", 
         "Pineal body": "Polymodal association"}

sois = list(sois_dict.keys())

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
#renaming for figure
short_nuclei = ["VA-L", "VM", "VPL", "VPM", "Post. Triangle", "PP", "Med. Geniculate", "dLGN", "LP", "Post. Complex", 
       "SGN", "Eth", "REth", "AV", "AM", "AD", "IAM", "IAD", "LD", "IMD",
       "MD", "Submedial", "PR", "Paraventricular", "PT", "Reuniens", "Xi", "RH", "CM", "PCN",
       "CL", "Parafascicular", "PIL", "RTN", "IntG", "vLGN", "SubG", "Med. Habenula", "Lat. Habenula", "PIN"]

sois_sort = np.array(short_nuclei)[order][:20]

#color palette based on nuclei type
cat = np.array(list(sois_dict.values())[1:])[order][:20] #removes thalamus soi
pal = [sns.color_palette("bright")[::-1][1] if n == "Sensory-motor" else sns.color_palette("bright")[::-1][0] for n in cat]

#boxplots of percent counts
plt.figure(figsize = (5,7))
df = pd.DataFrame(pcounts)
df.columns = short_nuclei
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort, palette = pal)
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={"facecolor":"None"}, order = sois_sort, palette = pal)
plt.xlabel("% of total thalamic cells")
plt.ylabel("Thalamic nuclei")

#make key
gold_patch = mpl.patches.Patch(color=sns.color_palette("bright")[::-1][1], label="Sensory-motor")
blue_patch = mpl.patches.Patch(color=sns.color_palette("bright")[::-1][0], label="Polymodal association")

plt.legend(title = "Thalamus nucleus type", 
           handles=[gold_patch, blue_patch], bbox_to_anchor=(.8, .9), loc=2, borderaxespad=0., frameon=False)
#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)

plt.savefig(os.path.join(dst, "thal_pcounts_boxplots.pdf"), bbox_inches = "tight")


#%%
#boxplots of density counts
order = np.argsort(np.median(density, axis = 0))[::-1]
#renaming for figure
short_nuclei = ["VA-L", "VM", "VPL", "VPM", "Post. Triangle", "PP", "Med. Geniculate", "dLGN", "LP", "Post. Complex", 
       "SGN", "Eth", "REth", "AV", "AM", "AD", "IAM", "IAD", "LD", "IMD",
       "MD", "Submedial", "PR", "Paraventricular", "PT", "Reuniens", "Xi", "RH", "CM", "PCN",
       "CL", "Parafascicular", "PIL", "RTN", "IntG", "vLGN", "SubG", "Med. Habenula", "Lat. Habenula", "PIN"]

sois_sort = np.array(short_nuclei)[order][:20]

#color palette based on nuclei type
cat = np.array(list(sois_dict.values())[1:])[order][:20] #removes thalamus soi
pal = [sns.color_palette("bright")[::-1][1] if n == "Sensory-motor" else sns.color_palette("bright")[::-1][0] for n in cat]

plt.figure(figsize = (5,7))
df = pd.DataFrame(density)
df.columns = short_nuclei
g = sns.stripplot(data = df,  color = "dimgrey", orient = "h", order = sois_sort, palette = pal)
sns.boxplot(data = df, orient = "h", showfliers=False, showcaps=False, 
            boxprops={"facecolor":"None"}, order = sois_sort, palette = pal)
plt.xlabel("Cells/$mm^3$")
plt.ylabel("Thalamic nuclei")

#make key
gold_patch = mpl.patches.Patch(color=sns.color_palette("bright")[::-1][1], label="Sensory-motor")
blue_patch = mpl.patches.Patch(color=sns.color_palette("bright")[::-1][0], label="Polymodal association")

plt.legend(title = "Thalamus nucleus type", 
           handles=[gold_patch, blue_patch], bbox_to_anchor=(.65, .7), loc=2, borderaxespad=0., frameon=False)
#hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)

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