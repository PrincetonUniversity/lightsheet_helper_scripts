#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:22:36 2019

@author: wanglab
"""

import pickle, numpy as np, matplotlib.pyplot as plt, os, matplotlib.ticker as ticker, pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from scipy import asarray as ar
from scipy.stats import chisquare

"""
You could try to use something like scipy.optimize.curve_fit to fit a 3D Gaussian, 
then calculate a chi^2 to assess the goodness of fit. You could see if there is a threshold in chi^2 
above which most of the edge cells lie and below which the human annotated cells lie.
What might be easier is seeing if there is a threshold in the width of the Gaussian that you fit.
"""
 
src = "/jukebox/wang/zahra/kelly_cell_detection_analysis"
ecells = os.path.join(src, "edge_cells.p")

rcells = os.path.join(src, "real_cells.p")

rcells = pickle.load(open(rcells, "rb"), encoding = "latin1")
ecells = pickle.load(open(ecells, "rb"), encoding = "latin1")

#real cells
#do for all x,y,z
yprof = []; cell_id = []; ints = []; df = pd.DataFrame()
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        cell_id.append((k,j))
        ints.append(cdct[j]["intensity"])
        yprof.append(cdct[j]["yprofile"])
#make np array
yprof = np.asarray(yprof)
#normalize
norm_yprof = []
for yy in yprof:    
    norm_yprof.append(np.asarray([(xx-np.nanmin(yy))/(np.nanmax(yy)-np.nanmin(yy)) for xx in yy]))    
norm_yprof = np.asarray(norm_yprof)

df["cell_id"] = cell_id

xprof = []; cell_id = []
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        xprof.append(cdct[j]["xprofile"])
#make np array
xprof = np.asarray(xprof)
#normalize
norm_xprof = []
for yy in xprof:    
    norm_xprof.append(np.asarray([(xx-np.nanmin(yy))/(np.nanmax(yy)-np.nanmin(yy)) for xx in yy]))
norm_xprof = np.asarray(norm_xprof)

zprof = []; cell_id = []; df_z = pd.DataFrame()
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        zprof.append(cdct[j]["zprofile"])
#make np array
zprof = np.asarray(zprof)
#normalize
norm_zprof = []
for yy in zprof:    
    norm_zprof.append(np.asarray([(xx-np.nanmin(yy))/(np.nanmax(yy)-np.nanmin(yy)) for xx in yy]))
norm_zprof = np.asarray(norm_zprof)

#do the same for edge cells
df_e = pd.DataFrame()
#do for all x,y,z
yprof_e = []; cell_id = []; ints_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        cell_id.append((k,j))
        ints_e.append(cdct[j]["intensity"])
        yprof_e.append(cdct[j]["yprofile"])
#make np array
yprof_e = np.asarray(yprof_e)
df_e["cell_id"] = cell_id
#normalize
norm_yprof_e = []
for yy in yprof_e:    
    norm_yprof_e.append(np.asarray([(xx-np.nanmin(yy))/(np.nanmax(yy)-np.nanmin(yy)) for xx in yy]))    
norm_yprof_e = np.asarray(norm_yprof_e)

xprof_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        xprof_e.append(cdct[j]["xprofile"])
#make np array
xprof_e = np.asarray(xprof_e)
#normalize
norm_xprof_e = []
for yy in xprof_e:    
    norm_xprof_e.append(np.asarray([(xx-np.nanmin(yy))/(np.nanmax(yy)-np.nanmin(yy)) for xx in yy]))
norm_xprof_e = np.asarray(norm_xprof_e)

zprof_e = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        zprof_e.append(cdct[j]["zprofile"])
#make np array
zprof_e = np.asarray(zprof_e)
#normalize
norm_zprof_e = []
for yy in zprof_e:    
    norm_zprof_e.append(np.asarray([(xx-np.nanmin(yy))/(np.nanmax(yy)-np.nanmin(yy)) for xx in yy]))
norm_zprof_e = np.asarray(norm_zprof_e)

#take mean
norm_z_mean = np.mean(norm_zprof, axis = 0)    
norm_y_mean = np.mean(norm_yprof, axis = 0)    
norm_x_mean = np.mean(norm_xprof, axis = 0)        


def get_chisq_pvals(profiles, norm, px=3, vis=False, cutoff=0.99):
    """ calculate chi square between mean cell profile and each cell profile """
    pvals = []; chistats = []
    for i in range(len(profiles)):
        cell = profiles[i][10-px:10+px+1]
        try:
            chistat, pval = chisquare(f_obs=cell, f_exp=norm[10-px:10+px+1])
            pvals.append(pval); chistats.append(chistat)
            if vis:
                if pval < cutoff:
                    plt.plot(norm[10-px:10+px+1], color = "blue")
                    plt.plot(cell, color = "blue", linestyle = "--")
        except:
            pval = -1
    pvals = np.asarray(pvals); chistats = np.asarray(chistats)
    
    return chistats, pvals

def get_cell_stats(profiles):
    diffs = []; mus = []; sigmas = []
    for i in range(len(profiles)):
        print(i)
        diffs.append(np.nanmax(profiles[i])-np.nanmin(profiles[i]))
        x = ar(range(len(profiles[i])))
        y = profiles[i]
        # weighted arithmetic mean (corrected - check the section below)
        mean = np.nansum(x * y) / np.nansum(y)
        sigma = np.sqrt(np.nansum(y * (x - mean)**2) / np.nansum(y))
        
        def Gauss(x, a, x0, sigma):
            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
        
        try:
            popt,pcov = curve_fit(Gauss, x, y, p0=[np.nanmax(y), mean, sigma])
            mx, mu, sigma = popt
        except Exception as e:
            print(e)
            mu, sigma = np.nan, np.nan
        mus.append(mu); sigmas.append(sigma)
        
    diffs = np.asarray(diffs); mus = np.asarray(mus); sigmas = np.asarray(sigmas)
    return diffs, mus, sigmas

#GET FEATURES
#uses normalized
#chisq fits for sliding window around cell
chistatsx_px3, pvalsx_px3 = get_chisq_pvals(norm_xprof, norm_x_mean, px = 3)
chistatsx_px5, pvalsx_px5 = get_chisq_pvals(norm_xprof, norm_x_mean, px = 5)
chistatsx_px10, pvalsx_px10 = get_chisq_pvals(norm_xprof, norm_x_mean, px = 10)

chistatsy_px3, pvalsy_px3 = get_chisq_pvals(norm_yprof, norm_y_mean, px = 3)
chistatsy_px5, pvalsy_px5 = get_chisq_pvals(norm_yprof, norm_y_mean, px = 5)
chistatsy_px10, pvalsy_px10 = get_chisq_pvals(norm_yprof, norm_y_mean, px = 10)

chistatsz_px3, pvalsz_px3 = get_chisq_pvals(norm_zprof, norm_z_mean, px = 3)
chistatsz_px5, pvalsz_px5 = get_chisq_pvals(norm_zprof, norm_z_mean, px = 5)
chistatsz_px10, pvalsz_px10 = get_chisq_pvals(norm_zprof, norm_z_mean, px = 10)

chistatsx_e_px3, pvalsx_e_px3 = get_chisq_pvals(norm_xprof_e, norm_x_mean, px = 3)
chistatsx_e_px5, pvalsx_e_px5 = get_chisq_pvals(norm_xprof_e, norm_x_mean, px = 5)
chistatsx_e_px10, pvalsx_e_px10 = get_chisq_pvals(norm_xprof_e, norm_x_mean, px = 10)

chistatsy_e_px3, pvalsy_e_px3 = get_chisq_pvals(norm_yprof_e, norm_y_mean, px = 3)
chistatsy_e_px5, pvalsy_e_px5 = get_chisq_pvals(norm_yprof_e, norm_y_mean, px = 5)
chistatsy_e_px10, pvalsy_e_px10 = get_chisq_pvals(norm_yprof_e, norm_y_mean, px = 10)

chistatsz_e_px3, pvalsz_e_px3 = get_chisq_pvals(norm_zprof_e, norm_z_mean, px = 3)
chistatsz_e_px5, pvalsz_e_px5 = get_chisq_pvals(norm_zprof_e, norm_z_mean, px = 5)
chistatsz_e_px10, pvalsz_e_px10 = get_chisq_pvals(norm_zprof_e, norm_z_mean, px = 10)

diffsx, mu_x, sigma_x = get_cell_stats(xprof)
diffsy, mu_y, sigma_y = get_cell_stats(yprof)
diffsz, mu_z, sigma_z = get_cell_stats(zprof)

diffsx_e, mu_xe, sigma_xe = get_cell_stats(xprof_e)
diffsy_e, mu_ye, sigma_ye = get_cell_stats(yprof_e)
diffsz_e, mu_ze, sigma_ze = get_cell_stats(zprof_e)

#fill dataframe with features
df["x_chisq_stat_px3"] = chistatsx_px3; df["y_chisq_stat_px3"] = chistatsy_px3; df["z_chisq_stat_px3"] = chistatsz_px3
df["x_chisq_pvals_px3"] = pvalsx_px3; df["y_chisq_pvals_px3"] = pvalsy_px3; df["z_chisq_pvals_px3"] = pvalsz_px3
#df["x_chisq_stat_px5"] = chistatsx_px5; df["y_chisq_stat_px5"] = chistatsy_px5; df["z_chisq_stat_px5"] = chistatsz_px5
#df["x_chisq_pvals_px5"] = pvalsx_px5; df["y_chisq_pvals_px5"] = pvalsy_px5; df["z_chisq_pvals_px5"] = pvalsz_px5
#df["x_chisq_stat_px10"] = chistatsx_px10; df["y_chisq_stat_px10"] = chistatsy_px10; df["z_chisq_stat_px10"] = chistatsz_px10
#df["x_chisq_pvals_px10"] = pvalsx_px10; df["y_chisq_pvals_px10"] = pvalsy_px10; df["z_chisq_pvals_px10"] = pvalsz_px10

df["x_diff_minima"] = diffsx; df["y_diff_minima"] = diffsy; df["z_diff_minima"] = diffsz
df["x_mean_guass"] = mu_x; df["y_mean_guass"] = mu_y; df["z_mean_guass"] = mu_z
df["x_sigma_guass"] = sigma_x; df["y_sigma_guass"] = sigma_y; df["z_sigma_guass"] = sigma_z

df_e["x_chisq_stat_px3"] = chistatsx_e_px3; df_e["y_chisq_stat_px3"] = chistatsy_e_px3; df_e["z_chisq_stat_px3"] = chistatsz_e_px3
df_e["x_chisq_pvals_px3"] = pvalsx_e_px3; df_e["y_chisq_pvals_px3"] = pvalsy_e_px3; df_e["z_chisq_pvals_px3"] = pvalsz_e_px3
#df_e["x_chisq_stat_px5"] = chistatsx_e_px5; df_e["y_chisq_stat_px5"] = chistatsy_e_px5; df_e["z_chisq_stat_px5"] = chistatsz_e_px5
#df_e["x_chisq_pvals_px5"] = pvalsx_e_px5; df_e["y_chisq_pvals_px5"] = pvalsy_e_px5; df_e["z_chisq_pvals_px5"] = pvalsz_e_px5
#df_e["x_chisq_stat_px10"] = chistatsx_e_px10; df_e["y_chisq_stat_px10"] = chistatsy_e_px10; df_e["z_chisq_stat_px10"] = chistatsz_e_px10
#df_e["x_chisq_pvals_px10"] = pvalsx_e_px10; df_e["y_chisq_pvals_px10"] = pvalsy_e_px10; df_e["z_chisq_pvals_px10"] = pvalsz_e_px10

df_e["x_diff_minima"] = diffsx_e; df_e["y_diff_minima"] = diffsy_e; df_e["z_diff_minima"] = diffsz_e
df_e["x_mean_guass"] = mu_xe; df_e["y_mean_guass"] = mu_ye; df_e["z_mean_guass"] = mu_ze
df_e["x_sigma_guass"] = sigma_xe; df_e["y_sigma_guass"] = sigma_ye; df_e["z_sigma_guass"] = sigma_ze

df.to_csv(os.path.join(src, "real_cell_stats.csv"), index = None)
df_e.to_csv(os.path.join(src, "edge_cell_stats.csv"), index = None)

