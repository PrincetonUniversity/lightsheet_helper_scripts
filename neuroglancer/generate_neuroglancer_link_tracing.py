#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 12:18:02 2019

@author: wanglab
"""

import neuroglancer
import cloudvolume
import pandas as pd

neuroglancer.set_static_content_source(url='https://neuromancer-seung-import.appspot.com')

# This volume handle can be used to notify the viewer that the data has changed.
viewer = neuroglancer.Viewer()

# Load in a segmentation layer (e.g. the raw-space allen atlas) that is being hosted with cloudvolume
with viewer.txn() as s:
    s.layers['20161205_tp_bl6_lob45_1000r_01'] = neuroglancer.ImageLayer(source='precomputed://http://localhost:1337',
    shader = '''
    void main() {
  float v = toNormalized(getDataValue(0)) * 20.0;
  emitRGBA(vec4(v, 0.0, 0.0, v));
}
'''                                                 
    )
print(viewer) # Click the link that is generated to see your volume displayed in Neuroglancer

with viewer.txn() as s:
    s.layers['rawatlas'] = neuroglancer.SegmentationLayer(source='precomputed://http://localhost:1338'
    )
    
# First you need to create a dictionary mapping the atlas id to the region name. 
# Here's my example

csv_file = '/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx'

df = pd.read_excel(csv_file)

atlas_df = df.copy()

atlas_dict = {}
for ind, row in atlas_df.iterrows():
    atlas_dict[int(row['id'])] = row['name']
    
# Here is the actual code for making the key binding. Copy and paste this
# into your jupyter notebook where you have already made the viewer() object

num_actions = 0
def my_action3(s):
    global num_actions
    num_actions += 1
    with viewer.config_state.txn() as st:
        region_id = s.selected_values['rawatlas']
        region_name = atlas_dict.get(region_id)
        st.status_messages['hello'] = ('%i:%s' %
                                    (region_id, region_name))

    print('Got my-action %i' % num_actions)
#     print('  Mouse position: %s' % (s.mouse_voxel_coordinates,))
    print('  Layer selected values: %s' % (s.selected_values,))
viewer.actions.add('my-action3', my_action3)

with viewer.config_state.txn() as s:
    s.input_event_bindings.viewer['keyp'] = 'my-action3'
    s.status_messages['hello'] = 'Welcome to the segment labeling example. Press "p" to see region name under cursor.'

ss = neuroglancer.ScreenshotSaver(viewer, '/home/wanglab/Documents/neuroglancer/screenshots')
ss.capture()   