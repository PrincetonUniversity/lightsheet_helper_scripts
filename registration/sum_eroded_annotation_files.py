# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:42:52 2020

@author: Zahra
"""

import os, tifffile, numpy as np

src = "/jukebox/wang/zahra/registration_error_pma/eroded_atlases/combined"
dst = os.path.join(src, "combined")

#make fld
if not os.path.exists(dst): os.mkdir(dst)

jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

lst = [xx for xx in os.listdir(src) if "tif" in xx]
lst.sort()

lst_job = lst[(jobid*5)-5:jobid*5]

print(lst_job)    

summed_ann = np.asarray([tifffile.imread(os.path.join(src, xx)) for xx in lst_job]).sum(axis = 0)

tifffile.imsave(os.path.join(dst, "annotation_pma_2018_20um_sagittal_erode_80um.tif"),
                summed_ann.astype("float32"))