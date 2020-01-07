#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:14:37 2019

@author: wanglab
"""

import os, subprocess as sp, tifffile, numpy as np, shutil, matplotlib.pyplot as plt, matplotlib as mpl
from collections import Counter

from tools.analysis.analyze_injection_inverse_transform import pool_injections_inversetransform
from tools.utils.io import makedir, load_kwargs, listdirfull
from tools.imageprocessing.orientation import fix_orientation
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


def run_transformix(invol, outpth, transformfile):
    
    #run transformix        
    sp.call(["transformix", "-in", invol, "-out", outpth, "-tp", transformfile])
    print("\n   Transformix File Generated: {}".format(outpth))
    
    return outpth

#here, im using the injection detection code for segmenting out the injection site from the resampled volume - more accurate
#later i can use that segmented volume in the downsampled space and convert it to atlas space

if __name__ == "__main__":

    src = "/jukebox/wang/pisano/tracing_output/retro_4x"
    
    brains = ["20180205_jg_bl6f_prv_01", "20180205_jg_bl6f_prv_02", "20180205_jg_bl6f_prv_03", "20180205_jg_bl6f_prv_04", 
          "20180215_jg_bl6f_prv_05", "20180215_jg_bl6f_prv_06",
       "20180215_jg_bl6f_prv_09", "20180305_jg_bl6f_prv_12", "20180305_jg_bl6f_prv_13","20180306_jg_bl6f_prv_14", 
       "20180305_jg_bl6f_prv_15", "20180306_jg_bl6f_prv_16", "20180312_jg_bl6f_prv_17", "20180326_jg_bl6f_prv_37",
       "20180313_jg_bl6f_prv_21", "20180313_jg_bl6f_prv_23", "20180313_jg_bl6f_prv_24", "20180305_jg_bl6f_prv_11", "20180313_jg_bl6f_prv_25",
       "20180322_jg_bl6f_prv_27", "20180322_jg_bl6f_prv_28", "20180323_jg_bl6f_prv_30", "20180326_jg_bl6f_prv_33", 
       "20180326_jg_bl6f_prv_34", "20180326_jg_bl6f_prv_35"]

    
    #run
    inputlist = [os.path.join(src, xx) for xx in brains]
    
    dct = {"inputlist": inputlist,
      "channel": "01",
      "channel_type": "injch",
      "filter_kernel": (3,3,3), 
      "threshold": 4, 
      "num_sites_to_keep": 1,
      "injectionscale": 45000, 
      "imagescale": 2,
      "reorientation": ("2","1","0"),
      "crop": "[:,:,:]",
      "dst": "/jukebox/wang/zahra/tracing_projects/prv/prv_injection_sites",
      "save_individual": False, 
      "save_tif": True,
      "colormap": "plasma", 
      "atlas": "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif",
      "annotation": "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif",
      "id_table": "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx"
    }
    
#    #only get brains for which the inj segmentation was successful
#    inj_brains = [os.path.join(src, xx[:-4]) for xx in os.listdir(dct["dst"]) if "tif" in xx]
#    
#    #even though we have the voxel counts in the csv file, i would still prefer to have the registered volumes just in case
#    #that is how we did the segmentation for the h129 anyways
#    #btw this will take long, and i don't recommend parallelizing bc of transformix
#    for fld in inj_brains:
#        
#        kwargs = load_kwargs(fld)
#        vol = [xx for xx in kwargs["volumes"] if xx.ch_type == "injch"][0] 
#        svlc = os.path.join(fld, "elastix")
#        #find transform file
#        sig2reg_fld = os.path.join(svlc, os.path.basename(vol.downsized_vol))+"/sig_to_reg"
#        transformfile = os.path.join(sig2reg_fld, "regtoatlas_TransformParameters.1.txt")
#        
#        #change the output image type bc otherwise it iterpolates too much and looks weird
#        with open(transformfile, "r") as file:
#            filedata = file.read()
#        # Replace the target string
#        filedata = filedata.replace('(ResultImagePixelType "short")', '(ResultImagePixelType "float")')
#        # Write the file out again
#        with open(transformfile, "w") as file:
#          file.write(filedata)
#    
#        invol = os.path.join(dct["dst"], os.path.basename(fld)+".tif")
#        outpth = os.path.join(dct["dst"], os.path.basename(fld)); makedir(outpth)
#        outpth = run_transformix(invol, outpth, transformfile)
#        
#        #fix the negative #'s around the site, and overwrite the tif onto the resampled version
#        #this was checked and is consistent w the transform
#        img = tifffile.imread(outpth+"/result.tif")
#        img[img < 0] = 0
#        assert np.sum(img < 0) == 0 #make sure there are no negative values
#        tifffile.imsave(invol, img.astype("uint16"))
#        shutil.rmtree(outpth) #delete original transformed file
        
    #inspect injection sites for the brains i currently have
    imgs = [os.path.join(dct["dst"], xx+".tif") for xx in brains]; imgs.sort()
    sites = np.array([fix_orientation(tifffile.imread(xx)[:, 450:, :], dct["reorientation"]) for xx in imgs]) #the y-axis cutoff for visualization
    #check
    for i,site in enumerate(sites):
        tifffile.imsave("/home/wanglab/Desktop/{}.tif".format(brains[i]), np.max(site, axis = 0))
    atl = fix_orientation(tifffile.imread(dct["atlas"])[:, 450:, :], dct["reorientation"])
    
    #make counter
    nonzeros = [list(zip(*np.nonzero(site))) for site in sites] #<-for pooled image
    
    #condense nonzero pixels
    nzs = [str(x) for xx in nonzeros for x in xx] #this list has duplicates if two brains had the same voxel w label
    c = Counter(nzs)
    array = np.zeros(atl.shape)
    print("Collecting nonzero pixels for pooled image...")
    tick = 0
    #generating pooled array where voxel value = total number of brains with that voxel as positive
    for k,v in c.items():
        k = [int(xx) for xx in k.replace("(","").replace(")","").split(",")]
        array[k[0], k[1], k[2]] = int(v)
        tick+=1
        if tick % 50000 == 0: print("   {}".format(tick))
        
    my_cmap = eval("plt.cm.{}(np.arange(plt.cm.RdBu.N))".format("viridis"))
    my_cmap[:1,:4] = 0.0  
    my_cmap = mpl.colors.ListedColormap(my_cmap)
    my_cmap.set_under("w")
    plt.figure()
    plt.imshow(np.max(atl, axis=0), cmap="gray")
    plt.imshow(np.max(array, axis=0), alpha=0.90, cmap=my_cmap); plt.colorbar(); plt.axis("off")
    
    plt.savefig(os.path.join("/home/wanglab/Desktop/heatmap_horz.pdf"), dpi = 300, transparent = True);
    plt.close()
#%%
    ann_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif")    
    #percentage covered by at least one injection plot (for supp figure)
    #annotation IDs of the cerebellum ONLY that are actually represented in annotation file
    ann = fix_orientation(ann_raw[:, 450:, :], dct["reorientation"])
    
    iids = {"Lingula (I)": 912,
            "Lobule II": 976,
            "Lobule III": 984,
            "Lobule IV-V": 1091,
            "Lobule VIa": 936,
            "Lobule VIb": 1134,
            "Lobule VII": 944,
            "Lobule VIII": 951,
            "Lobule IX": 957, #uvula IX
            "Lobule X": 968, #nodulus X
            "Simplex lobule": 1007, #simplex
            "Crus 1": 1056, #crus 1
            "Crus 2": 1064, #crus 2
            "Paramedian lobule": 1025, #paramedian lob
            "Copula pyramidis": 1033 #copula pyramidis
            }
    ak = np.asarray([k for k,v in iids.items()])
    
    atlas_rois = {}
    for nm, iid in iids.items():
        z,y,x = np.where(ann == iid) #find where structure is represented
        ann_blank = np.zeros_like(ann)
        ann_blank[z,y,x] = 1 #make a mask of structure in annotation space
        atlas_rois[nm] = ann_blank.astype(bool) #add mask of structure to dictionary
    
    expr_all_as_frac_of_lob = np.array([(array.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for nm, lob in atlas_rois.items()])    
    primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])


