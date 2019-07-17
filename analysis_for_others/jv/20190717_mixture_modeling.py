#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 12:53:45 2019

@author: wanglab
"""

from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from scipy.stats import multivariate_normal as mvn
from sklearn.preprocessing import StandardScaler

#read data
src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/pooled_analysis/60um_erosion_analysis/more_selected_structures/select_structures_percent_counts_for_plots.csv"

df = pd.read_csv(src)

#PCA

#make dataframe of percent counts per condition for signfiicant structures
dreadds = pd.DataFrame()
cnorev = pd.DataFrame()
cnonorev = pd.DataFrame()

structures = df.name.unique()

for soi in structures:
    #DREADDs
    dreadds[soi] = pd.Series(df[(df.name == soi)& (df.condition == "DREADDs")].percent.values)
    dreadds["condition"] = pd.Series(df[(df.name == soi) & (df.condition == "DREADDs")].condition.values)
    #CNO reversal
    cnorev[soi] = pd.Series(df[(df.name == soi) & (df.condition == "CNO_control_reversal")].percent.values)
    cnorev["condition"] = pd.Series(df[(df.name == soi) & (df.condition == "CNO_control_reversal")].condition.values)
    #CNO no reversal
    cnonorev[soi] = pd.Series(df[(df.name == soi) & (df.condition == "CNO_control_no_reversal")].percent.values)
    cnonorev["condition"] = pd.Series(df[(df.name == soi) & (df.condition == "CNO_control_no_reversal")].condition.values)
    
df_PCA = pd.concat([dreadds, cnorev, cnonorev])

#%%
# Separating out the features
arr = df_PCA.loc[:, structures].values
x = arr

# Separating out the target
y = df_PCA.loc[:,["condition"]].values

# Standardizing the features
x = StandardScaler().fit_transform(x)

#clean data
from sklearn.preprocessing import Imputer
imp = Imputer(missing_values="NaN", strategy="mean", axis=1)
cleaned_data = imp.fit_transform(x)

sklearn_pca = PCA(n_components = 2)

Y_sklearn = sklearn_pca.fit_transform(cleaned_data)

gmm = GaussianMixture(n_components=3, covariance_type='full').fit(Y_sklearn)
prediction_gmm = gmm.predict(Y_sklearn)
probs = gmm.predict_proba(Y_sklearn)

centers = np.zeros((3,2))
for i in range(3):
    density = mvn(cov=gmm.covariances_[i], mean=gmm.means_[i]).logpdf(Y_sklearn)
    centers[i, :] = Y_sklearn[np.argmax(density)]

plt.figure(figsize = (10,8))
plt.scatter(Y_sklearn[:, 0], Y_sklearn[:, 1], c = prediction_gmm, s = 200, cmap='viridis')
plt.scatter(centers[:, 0], centers[:, 1], c = 'black', s=300, alpha=0.6);