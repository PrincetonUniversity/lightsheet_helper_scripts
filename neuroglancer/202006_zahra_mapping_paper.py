#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 13:17:57 2020

@author: zahra

code to run neuroglancer after generating precomputed volumes
goal is to overlay raw data + atlas and make videos/screenshots of key structures
before running scripts, activate lightsheet env in each window that has neuroglancer and cloud volume installed
make sure you are connected to Princeton VPN and mounted on scratch/bucket
"""
###WINDOW 1###
#in the first ipython window run:
import neuroglancer 
neuroglancer.set_static_content_source(url="https://nglancer.pni.princeton.edu")

###WINDOW 2###
#in a new ipython window:
from cloudvolume import CloudVolume
layer_dir = "/jukebox/scratch/zmd/save/contra_ipsi_projection_studies_20191125/20170207_db_bl6_crii_1300r_02/647"
vol = CloudVolume(f"file://{layer_dir}")
vol.viewer(port=1338)

###WINDOW 1###
#go back to first window and run
neuroglancer.set_static_content_source(url="https://nglancer.pni.princeton.edu")
viewer = neuroglancer.Viewer()
with viewer.txn() as s:
    s.layers["20170207_db_bl6_crii_1300r_02"] = neuroglancer.ImageLayer(source="precomputed://http://localhost:1338"
    )
print(viewer)
#this should add the above volume to the neuroglancer window

###WINDOW 3###
#to add another layer (aka the atlas), in a new ipython window:
from cloudvolume import CloudVolume
layer_dir = "/jukebox/scratch/zmd/save/contra_ipsi_projection_studies_20191125/20170207_db_bl6_crii_1300r_02/atlas"
vol = CloudVolume(f"file://{layer_dir}")
vol.viewer(port=1337) #make sure this port is different from the first    

###WINDOW 1###
#go back to first window and run
with viewer.txn() as s:
    s.layers["20170207_db_bl6_crii_1300r_02_atlas"] = neuroglancer.SegmentationLayer(source="precomputed://http://localhost:1337"
    )
print(viewer)
#this should add the atlas volume to the neuroglancer window

###WINDOW 1###
#take screenshots
import os
svdst = "/jukebox/wang/zahra/neuroglancer/screenshots/20170207_db_bl6_crii_1300r_02/cb_brainstem"
#make sure these directories exist
if not os.path.exists(os.path.dirname(svdst)): os.mkdir(os.path.dirname(svdst)) #brain directory
if not os.path.exists(svdst): os.mkdir(svdst) #structure directory
for i in range(314,824):
    if i%10==0: print(i)
    with viewer.config_state.txn() as s:
        s.show_ui_controls = False
        s.show_panel_borders = False
    with viewer.txn() as s:
        s.voxel_coordinates = [2934,5692,i] #the xy coords here are from the neuroglancer window
        #(where the L center scale is located)
    #optionally limit window size
#    with viewer.config_state.txn() as s:
#        s.viewer_size = [1000,1000]
    ss = neuroglancer.ScreenshotSaver(viewer, svdst)
    ss.capture(index=i)

#after, return controls to neuroglancer browser
with viewer.config_state.txn() as s: 
    s.show_ui_controls = True 
    s.show_panel_borders = True 
