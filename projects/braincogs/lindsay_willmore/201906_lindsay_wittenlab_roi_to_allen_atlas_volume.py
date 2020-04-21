#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 12:12:29 2019

@author: wanglab
"""

import os, numpy as np, time, cv2
from skimage.external import tifffile
import matplotlib.pyplot as plt
import matplotlib.colors
from skimage.morphology import ball

os.chdir("/jukebox/wang/zahra/lightsheet_copy")
from tools.conv_net.input.read_roi import read_roi_zip
from tools.registration.transform import transformed_pnts_to_allen_helper_func
from tools.registration.register import count_structure_lister
from tools.utils.io import load_kwargs

if __name__ == "__main__":
    
    src = "/jukebox/wang/willmore/lightsheet/20190510_fiber_placement"
    
    flds = [os.path.join(src, xx) for xx in os.listdir(src) if os.path.isdir(os.path.join(src, xx))]
    
    for fld in flds:
        start = time.time()
        
        brain = os.path.basename(fld)
        #load in ROIS - clicked in SAGITTAL volume
        roi_pth = fld + "/{}_20190602_fiber_points_RoiSet.zip".format(brain)
        #check if path exists obv
        if os.path.exists(roi_pth):
            
            #have to load kwargs...
            kwargs = load_kwargs(fld)
        
            #get rois
            zyx_rois_sag = np.asarray([[int(yy) for yy in xx.replace(".roi", "").split("-")] for xx in read_roi_zip(roi_pth, include_roi_name=True)])
            #slice into diff orientations
            zyx_rois_cor = np.asarray([[xx[1], xx[2], xx[0]] for xx in zyx_rois_sag])
            zyx_rois_hor = np.asarray([[xx[2], xx[1], xx[0]] for xx in zyx_rois_sag])
            
            #make destination path
            dst = os.path.join(fld, "points_merged_to_atlas")
            if not os.path.exists(dst): os.mkdir(dst)
            
            #export coordinates
            if os.path.exists(os.path.join(dst, "{}_allen_coordinates.txt".format(brain))): os.remove(os.path.join(dst, "{}_allen_coordinates.txt".format(brain)))
            with open(os.path.join(dst, "{}_allen_coordinates.txt".format(brain)), "a") as txt:
                txt.write("Allen Atlas CCF coordinates (zyx) in the horizontal orientation:\n%s\n" % zyx_rois_hor)
                txt.write("\nAllen Atlas CCF coordinates (zyx) in the saggital orientation:\n%s\n" % zyx_rois_sag)
                txt.write("\nAllen Atlas CCF coordinates (zyx) in the coronal orientation:\n%s" % zyx_rois_cor)
            
            #MERGED IMAGES TO ATLAS
            #atlas (sagittal)
            atl = tifffile.imread("/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif")
            atl_cnn = np.zeros_like(atl)
            
            #make a merged map for the volume as well   
            sag_site1 = zyx_rois_sag[0:5]; sag_site2 = zyx_rois_sag[5:];
            sag_site1_av = sag_site1.mean(axis = 0).astype(int); sag_site2_av = sag_site2.mean(axis = 0).astype(int)
            atl_cnn[sag_site1_av[0]-1:sag_site1_av[0]+1, sag_site1_av[1], sag_site1_av[2]] = 1
            atl_cnn[sag_site2_av[0]-1:sag_site2_av[0]+1, sag_site2_av[1], sag_site2_av[2]] = 1
            
            #apply dilation
            r = 3
            selem = ball(r)[int(r/2)]
            atl_cnn = atl_cnn.astype("uint8")
            atl_cnn = np.asarray([cv2.dilate(atl_cnn[i], selem, iterations = 1) for i in range(atl_cnn.shape[0])])
         
            merged_Allen = np.stack([atl, atl_cnn, np.zeros_like(atl)], -1)
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_sagittal.tif".format(brain)), merged_Allen)
            
            #MERGED IMAGES TO REGISTERED VOLUMES
            #get registered volumes for overlay - this will be a saggital image
            #see if there is the injection channel first
            img = [vol.ch_to_reg_to_atlas for vol in kwargs["volumes"] if vol.ch_type == "injch"]
            #if not just get the regular reg channel
            if len(img) == 0: img = [vol.reg_to_atlas_vol for vol in kwargs["volumes"] if vol.ch_type == "regch"]
            #just need the one in the list
            img = tifffile.imread(img[0])
            img_cnn = np.zeros_like(img)
            
            #make a merged map for the volume as well   
            img_cnn[sag_site1_av[0]-1:sag_site1_av[0]+1, sag_site1_av[1], sag_site1_av[2]] = 1
            img_cnn[sag_site2_av[0]-1:sag_site2_av[0]+1, sag_site2_av[1], sag_site2_av[2]] = 1
            
            #apply dilation
            r = 3
            selem = ball(r)[int(r/2)]
            img_cnn = img_cnn.astype("uint8")
            img_cnn = np.asarray([cv2.dilate(img_cnn[i], selem, iterations = 1) for i in range(img_cnn.shape[0])])
            
            #make the merged saggital stack
            merged_img_sag = np.stack([img, img_cnn, np.zeros_like(img)], -1)
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_registered_image_sagittal.tif".format(brain)), merged_img_sag)
            
            #MERGED CORONAL IMAGES TO REGISTERED VOLUMES AND ATLAS FOR MAXIP
            #reslice
            #both these images have the same dims
            coronal_Allen = np.transpose(atl, [1, 2, 0]) #make coronal sections
            coronal_img = np.transpose(img, [1, 2, 0])
            
            print(coronal_Allen.shape, coronal_img.shape)

            #init
            cor_cnn = np.zeros_like(coronal_img)
            
            #make new circles marking fiber site... use this to overlay both on atlas and image
            cor_site1 = zyx_rois_cor[0:5]; cor_site2 = zyx_rois_cor[5:];
            cor_site1_av = cor_site1.mean(axis = 0).astype(int); cor_site2_av = cor_site2.mean(axis = 0).astype(int)
            cor_cnn[cor_site1_av[0]-2:cor_site1_av[0]+2, cor_site1_av[1], cor_site1_av[2]] = 1
            cor_cnn[cor_site2_av[0]-2:cor_site2_av[0]+2, cor_site2_av[1], cor_site2_av[2]] = 1
            
            #apply dilation
            r = 3
            selem = ball(r)[int(r/2)]
            cor_cnn = cor_cnn.astype("uint8")
            cor_cnn = np.asarray([cv2.dilate(cor_cnn[i], selem, iterations = 1) for i in range(cor_cnn.shape[0])])
            
            #merge
            merged_Allen_coronal = np.stack([coronal_Allen, cor_cnn, np.zeros_like(coronal_Allen)], -1)
            merged_img_coronal = np.stack([coronal_img, cor_cnn, np.zeros_like(coronal_img)], -1)
            #merge both registered volume and atlas to see registration quality
            merged_img_Allen_coronal = np.stack([coronal_img, coronal_Allen, cor_cnn], -1)
            
            #next, make smaller sections to visualise site better
            z = np.nonzero(cor_cnn)[0]
            
            #find z range of sites
            site1 = z[0:(z.shape[0]/2)]; site2 = z[(z.shape[0]/2)+1:]
            zrange_site1 = range(min(site1)-1, max(site1)+2); zrange_site2 = range(min(site2)-1, max(site2)+2)
            
            #TO ATLAS
            #save out coronal sections - based on the fact that you click 5 points in each site in HORIZONTAL sections
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_and_registered_volume_coronal.tif".format(brain)), merged_img_Allen_coronal.astype("uint16"))
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_coronal.tif".format(brain)), merged_Allen_coronal)
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_coronal_site1_z{}_{}.tif".format(brain, min(zrange_site1), max(zrange_site1))), merged_Allen_coronal[zrange_site1])
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_Allen_coronal_site2_z{}_{}.tif".format(brain, min(zrange_site2), max(zrange_site2))), merged_Allen_coronal[zrange_site2])
            
            #doing a max projection, in case you just want to look at that
            maxip1 = np.max(merged_Allen_coronal[zrange_site1], 0)
            maxip2 = np.max(merged_Allen_coronal[zrange_site2], 0)
            
            alpha = 0.6 #determines transparency, don't need to alter
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #color makes cells pop
            plt.axis("off")
            plt.imshow(maxip1[...,0], "gist_yarg")
            plt.imshow(maxip1[...,1], cmap, alpha = alpha)
            plt.title("Points overlaid on Allen atlas", fontsize = "small")
            plt.savefig(os.path.join(dst, "{}_points_merged_to_Allen_coronal_site1_maxip_z{}_{}.pdf".format(brain, min(zrange_site1), max(zrange_site1))), dpi = 300)
            
            plt.imshow(maxip2[...,0], "gist_yarg")
            plt.imshow(maxip2[...,1], cmap, alpha = alpha)
            plt.title("Points overlaid on Allen atlas", fontsize = "small")

            plt.savefig(os.path.join(dst, "{}_points_merged_to_Allen_coronal_site2_maxip_z{}_{}.pdf".format(brain, min(zrange_site2), max(zrange_site2))), dpi = 300)
            
            #TO REGISTERED IMAGE
            #save out coronal sections - based on the fact that you click 5 points in each site in HORIZONTAL sections
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_registered_image_coronal.tif".format(brain)), merged_img_coronal)
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_registered_image_coronal_site1_z{}_{}.tif".format(brain, min(zrange_site1), max(zrange_site1))), merged_img_coronal[zrange_site1])
            tifffile.imsave(os.path.join(dst, "{}_points_merged_to_registered_image_coronal_site2_z{}_{}.tif".format(brain, min(zrange_site2), max(zrange_site2))), merged_img_coronal[zrange_site2])
            
            #doing a max projection, in case you just want to look at that
            maxip1 = np.max(merged_img_coronal[zrange_site1], 0)
            maxip2 = np.max(merged_img_coronal[zrange_site2], 0)
            
            alpha = 0.5 #determines transparency, don't need to alter
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"]) #color makes cells pop
            plt.axis("off")
            plt.imshow(maxip1[...,0], "gist_yarg")
            plt.imshow(maxip1[...,1], cmap, alpha = alpha)
            plt.title("Points overlaid on registered volume", fontsize = "small")
            plt.savefig(os.path.join(dst, "{}_points_merged_to_registered_image_coronal_site1_maxip_z{}_{}.pdf".format(brain, min(zrange_site1), max(zrange_site1))), dpi = 300)
            
            plt.imshow(maxip2[...,0], "gist_yarg")
            plt.imshow(maxip2[...,1], cmap, alpha = alpha)
            plt.title("Points overlaid on registered volume", fontsize = "small")
            plt.savefig(os.path.join(dst, "{}_points_merged_to_registered_image_coronal_site2_maxip_z{}_{}.pdf".format(brain, min(zrange_site2), max(zrange_site2))), dpi = 300)
            
            print("\n\ntook {} seconds to make merged maps for {}\n".format(time.time()-start, brain))
            
            #make allen structure LUT
            zyx_rois = zyx_rois_sag
            
            #convert to structure
            annotation_file = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_template_25_sagittal_forDVscans.tif"
            ann = tifffile.imread(annotation_file)
            points = transformed_pnts_to_allen_helper_func(list(zyx_rois), ann, order = "ZYX")    
            
            #make dataframe
            lut_path = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table.xlsx"
            df = count_structure_lister(lut_path, *points)
            df.to_excel(os.path.join(dst, "{}_allen_structures.xlsx".format(brain)))