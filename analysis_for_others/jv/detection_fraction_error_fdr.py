#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 16:59:33 2019

@author: wanglab
"""


import pandas as pd, numpy as np

data = pd.read_csv('/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/rank_anovas.csv')

nonsig = np.asarray([len([xx for xx in data['rank 0 anova pcounts pval'].dropna() if xx > cutoff and xx <= cutoff+0.05])
                    for cutoff in np.arange(0.05, 1, 0.05)])

#take average
nonsig_av = np.mean(nonsig)
numerator = nonsig_av
#calculate number of pvals in bin 1
denom = len([xx for xx in data['rank 0 anova pcounts pval'].values if xx <= 0.05])
fdr = numerator / denom
#find error
N=len(data['rank 0 anova pcounts pval'].dropna())
error = np.sqrt(fdr * (1-fdr)/N)

print(fdr)
print(error)
