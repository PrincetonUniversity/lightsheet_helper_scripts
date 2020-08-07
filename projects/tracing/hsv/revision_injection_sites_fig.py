#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 16:59:09 2020

@author: wanglab
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np,os

#brain 1: simplex injection
cellreg = "/jukebox/wang/pisano/tracing_output/antero_4x/20161205_tp_bl6_sim_750r_03/elastix/20161205_tp_bl6_sim_750r_03_4x_647_008na_1hfds_z7d5um_150msec_10polvp_resized_ch00/result.tif"
injreg = "/jukebox/wang/pisano/tracing_output/antero_4x/20161205_tp_bl6_sim_750r_03/elastix/20161205_tp_bl6_sim_750r_03_4x_488_555_051na_1hfds_z7d5um_50msec_10polvp_resized_ch01/result.tif"
dst = "/jukebox/wang/zahra/tracing_projects/mapping_paper/revision_images"
cell = tif.imread(cellreg)
inj = tif.imread(injreg)

cell[cell < 0] = 60000
inj[inj < 0] = 60000

cell = cell[:,400:,:]
inj = inj[:,400:,:]

plt.imshow(np.max(inj[300:350]*10, axis=0))
inj_cor = np.rot90(np.transpose(inj, [1, 0, 2]),axes=(2,1))
plt.imshow(np.max(inj_cor[50:100]*5, axis=0),cmap="gist_yarg")
plt.axis("off")
plt.savefig(os.path.join(dst, "20161205_tp_bl6_sim_750r_03_inj_reg_maxP_zpln450-500.jpg"), 
           bbox_inches="tight", dpi=300)

cell_cor = np.rot90(np.transpose(cell, [1, 0, 2]),axes=(2,1))
plt.imshow(np.max(cell_cor[50:100]*7, axis=0),cmap="gist_yarg")
plt.axis("off")
plt.savefig(os.path.join(dst, "20161205_tp_bl6_sim_750r_03_cell_reg_maxP_zpln450-500.jpg"), 
           bbox_inches="tight", dpi=300)

#%%
#brain 2: lob VI
injreg = "/home/wanglab/wang/pisano/tracing_output/antero_4x/20170410_tp_bl6_lob6a_ml_repro_01/elastix/20170410_tp_bl6_lob6a_ml_repro_01_488_555_051na_1hfds_z7d5um_50msec_10povlp_resized_ch01/result.tif"
cellreg = "/home/wanglab/wang/pisano/tracing_output/antero_4x/20170410_tp_bl6_lob6a_ml_repro_01/elastix/20170410_tp_bl6_lob6a_ml_repro_01_647_008na_1hfds_z7d5um_200msec_10povlp_resized_ch00/result.tif"
cell = tif.imread(cellreg)
inj = tif.imread(injreg)

cell[cell < 0] = 60000
inj[inj < 0] = 60000

cell = cell[:,400:,:]
inj = inj[:,400:,:]

inj_cor = np.rot90(np.transpose(inj, [1, 0, 2]),axes=(2,1))
plt.imshow(np.max(inj_cor[50:200]*3, axis=0),cmap="gist_yarg")
plt.axis("off")
plt.savefig(os.path.join(dst, "20170410_tp_bl6_lob6a_ml_repro_01_inj_reg_maxP_zpln450-600.jpg"), 
           bbox_inches="tight", dpi=300)

cell_cor = np.rot90(np.transpose(cell, [1, 0, 2]),axes=(2,1))
plt.imshow(np.max(cell_cor[100:200]*2.5, axis=0),cmap="gist_yarg")
plt.axis("off")
plt.savefig(os.path.join(dst, "20170410_tp_bl6_lob6a_ml_repro_01_cell_reg_maxP_zpln450-600.jpg"), 
           bbox_inches="tight", dpi=300)