#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:56:34 2019

@author: wanglab
"""

import os, cv2, matplotlib.pyplot as plt, matplotlib.patches as mpatches

pth = "/jukebox/wang/zahra/neuroglancer/screenshots/20170204_tp_bl6_cri_1750r_03/amygdala"
list_of_files = [os.path.join(pth, xx) for xx in os.listdir(pth) if "png" in xx]; list_of_files.sort()

#make pngs into video
dst = os.path.join(os.path.dirname(pth), '20170204_tp_bl6_cri_1750r_03_amyg.avi')

frame_array = []

fps = 10 #frame rate

for png in list_of_files:
    #reading each files
    img = cv2.imread(png)
    
    #add labels to video
    fig = plt.figure()
    
    #info from neuroglancer
    patch0 = mpatches.Patch(color="#C4791B", label="Basomedial amygdalar n., ant.")
    patch1 = mpatches.Patch(color="#6CFF50", label="Basomedial amygdalar n., post.")    
    patch2 = mpatches.Patch(color="#AF1398", label="Basolateral amygdalar n., ant.")
    patch3 = mpatches.Patch(color="#C4791B", label="Basolateral amygdalar n., post.")
    patch4 = mpatches.Patch(color="#1F4423", label="Central amygdalar n., cap.")
    patch5 = mpatches.Patch(color="#1F3FFF", label="Central amygdalar n., med.")
    patch6 = mpatches.Patch(color="#BBFF04", label="Central amygdalar n., lat.")
    patch7 = mpatches.Patch(color="#F8FF42", label="Anterior amygdala")
    patch8 = mpatches.Patch(color="#FF2B15", label="Medial amygdalar n.")
    patch9 = mpatches.Patch(color="#A14DFF", label="Lateral amygdalar n.")
    patch10 = mpatches.Patch(color="#E26CA8", label="Posterior amygdalar n.")
    patch11 = mpatches.Patch(color="#45B2FF", label="Intercalated amygdalar n.")
    patch12 = mpatches.Patch(color="#C4FF26", label="Cortical amygdala, ant.")
    patch13 = mpatches.Patch(color="#4AC1FF", label="Cortical amygdala, post., med.")
    patch14 = mpatches.Patch(color="#2659FF", label="Cortical amygdala, post., lat.")
                            
    plt.imshow(img)
    plt.legend(handles=[patch0, patch1, patch2, patch3, patch4, patch5, patch6, patch7, patch8,
                        patch9, patch10, patch11, patch12, patch13, patch14], loc=0, 
                        borderaxespad=0., fontsize=4, framealpha=0.4)
    plt.axis("off")    
    plt.savefig(png, dpi=300, bbox_inches = "tight") #overwrite file
    
    #read file again
    img = cv2.imread(png)
    height, width, layers = img.shape
    size = (width,height)
    
    #inserting the frames into an image array
    frame_array.append(img)
    plt.close()

out = cv2.VideoWriter(dst, cv2.VideoWriter_fourcc(*'DIVX'), fps, size)
    
for i in range(len(frame_array)):
    # writing to a image array
    out.write(frame_array[i])
out.release()