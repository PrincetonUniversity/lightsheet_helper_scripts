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

brainname = "20180322_jg_bl6f_prv_28"
###WINDOW 1###
#in the first ipython window run:
import neuroglancer 
neuroglancer.set_static_content_source(url="https://nglancer.pni.princeton.edu")
port=1337

###WINDOW 2###
#in a new ipython window:
from cloudvolume import CloudVolume
brainname = "20180322_jg_bl6f_prv_28"
port=1337
layer_dir = "/jukebox/scratch/zmd/save/contra_ipsi_projection_studies_20191125/%s/647" % brainname
vol = CloudVolume(f"file://{layer_dir}")
vol.viewer(port=port)

###WINDOW 1###
#go back to first window and run
neuroglancer.set_static_content_source(url="https://nglancer.pni.princeton.edu")
viewer = neuroglancer.Viewer()
with viewer.txn() as s:
    s.layers["%s" % brainname] = neuroglancer.ImageLayer(source="precomputed://http://localhost:%s" % port
    )
print(viewer)
#this should add the above volume to the neuroglancer window

###WINDOW 3###
#to add another layer (aka the atlas), in a new ipython window:
from cloudvolume import CloudVolume
brainname = "20180322_jg_bl6f_prv_28"
port=1337
layer_dir = "/jukebox/scratch/zmd/save/contra_ipsi_projection_studies_20191125/%s/atlas" % brainname
vol = CloudVolume(f"file://{layer_dir}")
vol.viewer(port=port+1) #make sure this port is different from the first    

###WINDOW 1###
#go back to first window and run
with viewer.txn() as s:
    s.layers["%s_atlas" % brainname] = neuroglancer.SegmentationLayer(source="precomputed://http://localhost:%s" % int(port+1)
    )
print(viewer)
#this should add the atlas volume to the neuroglancer window

###WINDOW 1###
#take screenshots
import os
svdst = "/jukebox/wang/zahra/neuroglancer/screenshots/%s/cb_brainstem" % brainname
#make sure these directories exist
if not os.path.exists(os.path.dirname(svdst)): os.mkdir(os.path.dirname(svdst)) #brain directory
if not os.path.exists(svdst): os.mkdir(svdst) #structure directory
for i in range(334,824):
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

#%%

import math, sys, matplotlib.pyplot as plt
# Define the functions for finding the color hex strings given the color seed that Neuroglancer sets
# and a segment id. These functions should not be modified.
def hash_function(state,value):
    """ Python implementation of hashCombine() function
    in src/neuroglancer/gpu_hash/hash_function.ts,
    a modified murmur hcolorash
    """
    k1 = 0xcc9e2d51
    k2 = 0x1b873593
    state = state & 0xffffffff
    value = (value * k1) & 0xffffffff
    value = ((value << 15) | value >> 17) & 0xffffffff
    value = (value * k2) & 0xffffffff
    state = (state ^ value) & 0xffffffff
    state = (( state << 13) | state >> 19) & 0xffffffff
    state = (( state * 5) + 0xe6546b64) & 0xffffffff
    return state

def hsv_to_rgb(h,s,v):
    """ Convert H,S,V values to RGB values.
    Python implementation of hsvToRgb in src/neuroglancer/util/colorspace.ts """
    h*=6
    hue_index = math.floor(h)
    remainder = h - hue_index
    val1 = v*(1-s)
    val2 = v*(1-(s*remainder))
    val3 = v*(1-(s*(1-remainder)))
    hue_remainder = hue_index % 6
    if hue_remainder == 0:
        return (v,val3,val1)
    elif hue_remainder == 1:
        return (val2,v,val1)
    elif hue_remainder == 2:
        return (val1,v,val3)
    elif hue_remainder == 3:
        return (val1,val2,v)
    elif hue_remainder == 4:
        return (val3,val1,v)
    elif hue_remainder == 5: 
        return (v,val1,val2)   

def pack_color(rgb_vec):
    """ Returns an integer formed
    by concatenating the channels of the input color vector.
    Python implementation of packColor in src/neuroglancer/util/color.ts
    """
    result = 0
    for i in range(len(rgb_vec)):
        result = ((result << 8) & 0xffffffff) + min(255,max(0,round(rgb_vec[i]*255)))
    return result

def hex_string_from_segment_id(color_seed,segment_id):
    """ Return the hex color string for a segment
    given a color seed and the segment id """
    segment_id = int(segment_id) # necessary since segment_id is 64 bit originally 
    result = hash_function(state=color_seed,value=segment_id)
    newvalue = segment_id >> 32
    result2 = hash_function(state=result,value=newvalue)
    c0 = (result2 & 0xFF) / 255.
    c1 = ((result2 >> 8) & 0xFF) / 255.;
    h = c0
    s =  0.5 + 0.5 * c1
    v = 1.0
    rgb=hsv_to_rgb(h,s,v)
    packed_color = pack_color(rgb_vec=rgb)
    hex_string = format(packed_color, 'x')
    """ Zero pad the hex string if less than 6 characeters """
    if len(hex_string) < 6:
        hex_string = '0'*(6-len(hex_string)) + hex_string
    hex_string = '#' + hex_string
    return hex_string

# Gets the list of active segments, the current color seed and any manual changes 
with viewer.txn() as s:
    seglayer = s.layers['atlas'].layer
    seg_dict = seglayer.to_json()
    try:
        color_seed = seg_dict['colorSeed']
    except:
        color_seed = 0
    try:
        active_segments=list(seg_dict['segments'])
    except:
        active_segments=[]
    try:
        manual_segment_dict = seg_dict['segmentColors']
    except:
        manual_segment_dict = {}
if active_segments == []:
    sys.exit("You do not have any active segments selected. Select some segments and try re-running this block")
    
# Fill a dictionary where
# keys will be segment id, values will be hex color string
# taking the manually set value over the original one if present
hex_str_dict = {}
for segment_id in active_segments:
    if segment_id in manual_segment_dict:
        hex_str_dict[segment_id] = manual_segment_dict[segment_id]
    else:
        hex_str_dict[segment_id] = hex_string_from_segment_id(color_seed,segment_id)
print("Colors of your segments are:")
print(hex_str_dict)
print()

sizes = [360/float(len(active_segments)) for x in active_segments]

fig1, ax1 = plt.subplots()
pie=ax1.pie(sizes,labels=hex_str_dict.keys(),labeldistance=1.0,
        textprops={'fontsize': 12}, startangle=90,colors=hex_str_dict.values())
title=fig1.suptitle('Verify segment colors:',fontsize=18)