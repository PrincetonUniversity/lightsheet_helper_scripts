#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:56:34 2019

@author: wanglab
"""

import os, cv2, matplotlib.pyplot as plt, matplotlib.patches as mpatches

if __name__ == "__main__":
    
    pth = "/jukebox/wang/zahra/neuroglancer/screenshots/20170115_tp_bl6_lob6b_ml_04/frontal_areas"
    dst = pth+"_w_labels"
    if not os.path.exists(dst): os.mkdir(dst)
    list_of_files = [os.path.join(pth, xx) for xx in os.listdir(pth) if "png" in xx]; list_of_files.sort()
    #make pngs into video
    vid_dst = os.path.join(os.path.dirname(pth), '20170115_tp_bl6_lob6b_ml_04_pfc.avi')
    
    frame_array = []
    
    fps = 10 #frame rate
    
    alpha = 0.4 #transparency
    
    for png in list_of_files:
        #reading each files
        img = plt.imread(png)
        #add labels to video
        plt.figure()
        #info from neuroglancer
        patch0 = mpatches.Patch(color="#4EB829", label="Anterior cingulate, ventral", alpha=alpha)
        patch1 = mpatches.Patch(color="#EDFF75", label="Anterior cingulate, dorsal", alpha=alpha)  
        patch2 = mpatches.Patch(color="#6036FF", label="Prelimbic", alpha=alpha)
        patch3 = mpatches.Patch(color="#1CDDEC", label="Orbital, medial", alpha=alpha)
        patch4 = mpatches.Patch(color="#C56835", label="Orbital, ventrolateral", alpha=alpha)
        patch5 = mpatches.Patch(color="#D38DED", label="Orbital, lateral", alpha=alpha)
        patch6 = mpatches.Patch(color="#DF274A", label="Infralimbic", alpha=alpha)
    #    patch7 = mpatches.Patch(color="#F8FF42", label="Anterior amygdala")
    #    patch8 = mpatches.Patch(color="#FF2B15", label="Medial amygdalar n.")
    #    patch9 = mpatches.Patch(color="#A14DFF", label="Lateral amygdalar n.")
    #    patch10 = mpatches.Patch(color="#E26CA8", label="Posterior amygdalar n.")
    #    patch11 = mpatches.Patch(color="#45B2FF", label="Intercalated amygdalar n.")
    #    patch12 = mpatches.Patch(color="#C4FF26", label="Cortical amygdala, ant.")
    #    patch13 = mpatches.Patch(color="#4AC1FF", label="Cortical amygdala, post., med.")
    #    patch14 = mpatches.Patch(color="#2659FF", label="Cortical amygdala, post., lat.")
        plt.imshow(img)
        plt.legend(handles=[patch0, patch1, patch2, patch3, patch4, patch5, patch6],
    #                        , patch7, patch8,
    #                        patch9, patch10, patch11, patch12, patch13, patch14], 
                            loc=4, 
                            borderaxespad=0., fontsize=5, framealpha=0.4)
        plt.axis("off")    
        plt.savefig(os.path.join(dst, os.path.basename(png)), dpi=300, bbox_inches = "tight", pad_inches = 0) #write file to new location
        plt.close()
        
        
    #now make labeled images into videos
    list_of_lbld_files = [os.path.join(dst, xx) for xx in os.listdir(dst) if "png" in xx]; list_of_lbld_files.sort()
    
    for png in list_of_lbld_files:
        #read file again
        img = cv2.imread(png)
        height, width, layers = img.shape
        size = (width,height)
        
        #inserting the frames into an image array
        frame_array.append(img)
    
    out = cv2.VideoWriter(vid_dst, cv2.VideoWriter_fourcc(*'DIVX'), fps, size)
        
    for i in range(len(frame_array)):
        # writing to a image array
        out.write(frame_array[i])
    out.release()