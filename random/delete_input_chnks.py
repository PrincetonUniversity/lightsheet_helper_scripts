#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 11:35:01 2018

@author: wanglab
"""

import os, shutil

scratch_dir = "/jukebox/scratch/zmd"

#check if brain has completed inference and reconstructed array exists
input_chnks_delete = [xx for xx in os.listdir(scratch_dir) if "output_chnks" in os.listdir(os.path.join(scratch_dir, xx))
                        and "input_chnks" in os.listdir(os.path.join(scratch_dir, xx))                        
                        and "reconstructed_array.npy" in os.listdir(os.path.join(scratch_dir, xx)) 
                        and len(os.listdir(os.path.join(os.path.join(scratch_dir, xx), "input_chnks"))) == 
                        len(os.listdir(os.path.join(os.path.join(scratch_dir, xx), "output_chnks")))]

for brain in input_chnks_delete:
    shutil.rmtree(os.path.join(os.path.join(scratch_dir, brain), "input_chnks"))
    print("cnn inference and reconstruction run successfully!\
          \ninput chunks deleted for: \n{}".format(os.path.join(scratch_dir, brain)))
    
    
#[zmd@spock-login python]$ python delete_input_chnks.py 
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20170115_tp_bl6_lob6b_500r_05
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20170116_tp_bl6_lob6b_lpv_07
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20170115_tp_bl6_lob6a_rpv_03
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20170308_tp_bl6f_lob7_2x_02
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20170308_tp_bl6_lob8_lat_05
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180327_jg40_bl6_sim_03
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180327_jg42_bl6_lob6a_05
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180410_jg49_bl6_lob45_02
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180409_jg46_bl6_lob6a_04
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180409_jg47_bl6_lob6a_05
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180409_jg44_bl6_lob6a_02
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180410_jg50_bl6_lob6b_03
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180409_jg45_bl6_lob6a_03
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180410_jg51_bl6_lob6b_04
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180409_jg43_bl6_lob6a_01
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180410_jg48_bl6_lob6a_01
#cnn inference and reconstruction run successfully!          
#input chunks deleted for: 
#/jukebox/scratch/zmd/20180410_jg52_bl6_lob7_05
