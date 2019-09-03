#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 12:20:34 2019

@author: wanglab
"""

import os, scipy.io as sio, numpy as np, tifffile as tif, pandas as pd, matplotlib.pyplot as plt, cv2, seaborn as sns
from skimage.morphology import ball
from scipy.stats import chisquare
from scipy import asarray as ar
from scipy.optimize import curve_fit
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer

#setup
img_pth = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos"
brainname = "171213_m37110_demonstrator_20171016_790_015na_1hfsds_z5um_1000msec_16-03-36"
src = "/jukebox/wang/zahra/kelly_cell_detection_analysis"
brain = os.path.join(img_pth, brainname)

dct = sio.loadmat(os.path.join(brain, "sliding_diff_peak_find_975percentile_20190227_format2.mat"))
arr = dct["cell_centers_orig_coord"]

def get_features_from_cell(cell, vol, w=10, px=3):
    
    #based on labeled data, obv figure out a way to import them in the future
    norm_z_mean = [0.20293142, 0.21529837, 0.23150973, 0.25075755, 0.2467855 ,
           0.28805012, 0.3697442 , 0.511728  , 0.68492967, 0.82598543,
           0.97048056, 0.82300526, 0.6743533 , 0.49407753, 0.34736645,
           0.24914072, 0.20730427, 0.19215529, 0.18126972, 0.17916188,
           0.1542593]
    norm_y_mean = [0.19051766, 0.198287  , 0.19966966, 0.20214507, 0.197877  ,
           0.1979725 , 0.2073713 , 0.23559561, 0.3275807 , 0.5976898 ,
           0.9503866 , 0.6140519 , 0.31939164, 0.23095436, 0.20971417,
           0.21110429, 0.21813215, 0.23384634, 0.24211957, 0.24829383,
           0.2459892 ]
    norm_x_mean = [0.17778772, 0.1984067 , 0.19856174, 0.19846928, 0.20103398,
           0.18993343, 0.2106994 , 0.2570913 , 0.37046623, 0.6107595 ,
           0.9397369 , 0.5924273 , 0.3130014 , 0.22176723, 0.20720865,
           0.22059521, 0.22673169, 0.23923236, 0.23971424, 0.24849065,
           0.25069913]

    z,y,x = cell; z=z-700; y=y-1; x=x-1 #offsets
    xprofile = vol[z,y,x-w:x+w+1] #cell centered at index = 3
    yprofile = vol[z,y-w:y+w+1,x] #cell centered at index = 3
    zprofile = vol[z-w:z+w+1,y,x] #cell centered at index = 3
    intensity = vol[z,y,x]
    chistatx, pvalx = chisquare(f_obs=[((xx-min(xprofile))/(max(xprofile)-min(xprofile))) for xx in xprofile][10-px:10+px+1],  #normalized
                                       f_exp=norm_x_mean[10-px:10+px+1])
    chistaty, pvaly = chisquare(f_obs=[((xx-min(yprofile))/(max(yprofile)-min(yprofile))) for xx in yprofile][10-px:10+px+1], 
                                       f_exp=norm_y_mean[10-px:10+px+1])
    chistatz, pvalz = chisquare(f_obs=[((xx-min(zprofile))/(max(zprofile)-min(zprofile))) for xx in zprofile][10-px:10+px+1], 
                                       f_exp=norm_z_mean[10-px:10+px+1])
    
    #difference bw minimas
    diffx, diffy, diffz = np.max(xprofile)-np.min(xprofile), np.max(yprofile)-np.min(yprofile), np.max(zprofile)-np.min(zprofile)
    
    def get_guassian_stats(profile):
        x = ar(range(len(profile)))
        y = profile
        mean = np.sum(x * y) / np.sum(y)
        sigma = np.sqrt(np.sum(y * (x - mean)**2) / np.sum(y))
            
        def Gauss(x, a, x0, sigma):
            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
        
        try:
            popt,pcov = curve_fit(Gauss, x, y, p0=[np.nanmax(y), mean, sigma])
            mx, mu, sigma = popt
        except Exception as e:
            print(e)
            mu, sigma = np.nan, np.nan
        
        return mu, sigma
   
    mux, sigmax = get_guassian_stats(xprofile)
    muy, sigmay = get_guassian_stats(yprofile)
    muz, sigmaz = get_guassian_stats(zprofile)
    
    return chistatx, chistaty, chistatz, pvalx, pvaly, pvalz, diffx, diffy, diffz, mux, muy, muz, sigmax, sigmay, sigmaz, intensity
    
#lets do some middle 30 planes
zrange = np.arange(700, 730)
imgs = [os.path.join(brain, xx) for yy in zrange for xx in os.listdir(brain) if "tif" in xx and "Z%04d"%yy in xx]; imgs.sort()
vol = np.asarray([tif.imread(img) for img in imgs]) #image volumes

cells = arr[(arr[...,0] >= 710) & (arr[...,0] < 720)] #cells detected in image volume

#extract all the features from these cells
features = np.zeros((len(cells), 16)) #16 features
for i,cell in enumerate(cells):
    print(i)
    features[i] = get_features_from_cell(cell, vol)
    
 
rcells = pd.read_csv(os.path.join(src, "real_cell_stats.csv"))
ecells = pd.read_csv(os.path.join(src, "edge_cell_stats.csv"))

#find variables
rc = rcells.copy().dropna()
params = [xx for xx in rc.columns if xx != "cell_id"]

#combine both cells into one dataset with a label colum
ecs = ecells.drop(columns=["cell_id"]).dropna()
rcs = rcells.drop(columns=["cell_id"]).dropna()
ecs["label"] = np.zeros(len(ecs))
rcs["label"] = np.ones(len(rcs))
allcells = pd.concat([ecs, rcs])

X = allcells[params] # Features
y = allcells.label # Target variable

X_train,X_test,y_train,y_test=train_test_split(X,y,test_size=0.25)
# instantiate the model (using the default parameters)
logreg = LogisticRegression()
# fit the model with data
logreg.fit(X_train,y_train)

#predict on real brain data
#impute missing values
imp = SimpleImputer(missing_values=np.nan, strategy='median')
imp = imp.fit(features)
features_imp = imp.transform(features)
y_pred = logreg.predict(features_imp)

#alternatively, drop guassian values?
#features_imp = features[...,:9]
y_pred = logreg.predict(features_imp)
#%%
#visualize results
classified_cells = cells[y_pred.astype(bool)]

prob = logreg.predict_proba(features_imp)
sns.boxplot(prob[...,1], orient = "h")
plt.xlabel("Classifier probability for real cells")
plt.ylabel("Cells (n = {})".format(len(prob)))
plt.axvline(x=0.5, color = "gray", linestyle = "--")
#%%
#map a cell map
cell_map = np.zeros_like(vol)
cell_map_c = np.zeros_like(vol)

#before classifier
for cell in cells:
    z,y,x = cell
    z=z-700; y=y-1; x=x-1 #adjusting offsets
    cell_map[z,y,x] = 6000

for cell in classified_cells:
    z,y,x = cell
    z=z-700; y=y-1; x=x-1 #adjusting offsets
    cell_map_c[z,y,x] = 6000

#apply x y dilation
r = 2
selem = ball(r)[int(r/2)]
cell_map = cell_map.astype("uint8"); cell_map_c = cell_map_c.astype("uint8")
cell_map = np.asarray([cv2.dilate(cell_map[i], selem, iterations = 1) for i in range(cell_map.shape[0])])
cell_map_c = np.asarray([cv2.dilate(cell_map_c[i], selem, iterations = 1) for i in range(cell_map_c.shape[0])])

merged = np.stack([vol, cell_map, cell_map_c], -1) #rgb image you can open up in fiji; volume = red; cells = green

tif.imsave(os.path.join(src, "DV_before_after_classfier.tif"), merged)