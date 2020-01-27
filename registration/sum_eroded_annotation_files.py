# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:42:52 2020

@author: Zahra
"""

import os, tifffile, numpy as np

src = '/jukebox/wang/zahra/kelly_cell_detection_analysis/eroded_atlases/combined'
#src = r'Z:\zahra\kelly_cell_detection_analysis\eroded_atlases'

#jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

lst = os.listdir(src)
lst.sort()

#lst_job = lst[(jobid*10)-10:jobid*10]

#print(lst_job)    

summed_ann = np.asarray([tifffile.imread(os.path.join(src, xx)) for xx in lst]).sum(axis = 0)

tifffile.imsave("/jukebox/wang/zahra/kelly_cell_detection_analysis/annotation_allen_2017_25um_sagittal_erode_80um.tif",
                summed_ann.astype("uint16"))