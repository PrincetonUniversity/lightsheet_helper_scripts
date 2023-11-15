#! /bin/env python

import os, sys

from cloudvolume import CloudVolume
from cloudvolume.lib import mkdir, touch
import shutil
import json
# for Princeton mouse atlas:
#src_seg_props_info = '/jukebox/LightSheetData/atlas/neuroglancer/atlas/princetonmouse/segment_properties/info'
# for allen atlas:
src_seg_props_info = '/jukebox/LightSheetTransfer/atlas/neuroglancer/atlas/allenatlas_2017/segment_properties/info'
# Set viz_dir to parent folder where you put your precomputed raw atlas layers
viz_dir="/jukebox/LightSheetData/lightserv_testing/neuroglancer/jess_cfos" # bucket dir, so layers will be in $BUCKET/$DATASET/$LAYER 

brains = [] # put list of brain ids here, e.g. "f37077_demonstrator"
if __name__ == "__main__":
    brains = [] # put list of animal ids in here
    for brain in brains:
            # Update the info file and save it
            print(f"Dataset: {dataset}, animal_id: {brain}")
            layer_dir = os.path.join(viz_dir,f'{brain}_raw_atlas')
            vol = CloudVolume(f'file://{layer_dir}')
            info_dict = vol.info
            info_dict['segment_properties'] = "segment_properties"
            info_filename = '/'.join(vol.info_cloudpath.split('/')[2:]) 
            with open(info_filename,'w') as outfile:
                json.dump(info_dict,outfile,sort_keys=True,indent=2)
            print(f"ammended info file to include 'segment_properties' key: {info_filename}")
            # copy over the segment_properties directory
            seg_props_dir = os.path.join(layer_dir,'segment_properties')
            mkdir(seg_props_dir)
            dest_seg_props_info = os.path.join(seg_props_dir,'info')
            shutil.copyfile(src_seg_props_info,dest_seg_props_info)
            print("copied over segment_properties info file")
