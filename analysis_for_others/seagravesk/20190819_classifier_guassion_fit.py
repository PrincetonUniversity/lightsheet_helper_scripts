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

dst = "/jukebox/wang/zahra/kelly_cell_detection_analysis"


#look at real cells first
#do for all x,y,z
yprof = []; cell_id = []; df = pd.DataFrame()
#make 1 array for all y profiles of cells
for k,v in rcells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        cell_id.append((k,j))
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
yprof_e = []; cell_id = []
#make 1 array for all y profiles of cells
for k,v in ecells.items():
    cdct = v
    for j,m in cdct.items():
        print(j)
        cell_id.append((k,j))
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

#uses normalized
chistatsx, pvalsx = get_chisq_pvals(norm_xprof, norm_x_mean)
chistatsy, pvalsy = get_chisq_pvals(norm_yprof, norm_y_mean)
chistatsz, pvalsz = get_chisq_pvals(norm_zprof, norm_z_mean)
chistatsx_e, pvalsx_e = get_chisq_pvals(norm_xprof_e, norm_x_mean)
chistatsy_e, pvalsy_e = get_chisq_pvals(norm_yprof_e, norm_y_mean)
chistatsz_e, pvalsz_e = get_chisq_pvals(norm_zprof_e, norm_z_mean)
diffsx, mu_x, sigma_x = get_cell_stats(xprof)
diffsy, mu_y, sigma_y = get_cell_stats(yprof)
diffsz, mu_z, sigma_z = get_cell_stats(zprof)
diffsx_e, mu_xe, sigma_xe = get_cell_stats(xprof_e)
diffsy_e, mu_ye, sigma_ye = get_cell_stats(yprof_e)
diffsz_e, mu_ze, sigma_ze = get_cell_stats(zprof_e)
df["x_chisq_stat"] = chistatsx; df["y_chisq_stat"] = chistatsy; df["z_chisq_stat"] = chistatsz
df["x_chisq_pvals"] = pvalsx; df["y_chisq_pvals"] = pvalsy; df["z_chisq_pvals"] = pvalsz
df["x_diff_minima"] = diffsx; df["y_diff_minima"] = diffsy; df["z_diff_minima"] = diffsz
df["x_mean_guass"] = mu_x; df["y_mean_guass"] = mu_y; df["z_mean_guass"] = mu_z
df["x_sigma_guass"] = sigma_x; df["y_sigma_guass"] = sigma_y; df["z_sigma_guass"] = sigma_z

df_e["x_chisq_stat"] = chistatsx_e; df_e["y_chisq_stat"] = chistatsy_e; df_e["z_chisq_stat"] = chistatsz_e
df_e["x_chisq_pvals"] = pvalsx_e; df_e["y_chisq_pvals"] = pvalsy_e; df_e["z_chisq_pvals"] = pvalsz_e
df_e["x_diff_minima"] = diffsx_e; df_e["y_diff_minima"] = diffsy_e; df_e["z_diff_minima"] = diffsz_e
df_e["x_mean_guass"] = mu_xe; df_e["y_mean_guass"] = mu_ye; df_e["z_mean_guass"] = mu_ze
df_e["x_sigma_guass"] = sigma_xe; df_e["y_sigma_guass"] = sigma_ye; df_e["z_sigma_guass"] = sigma_ze

df.to_csv(os.path.join(dst, "real_cell_stats.csv"), index = None)
df_e.to_csv(os.path.join(dst, "edge_cell_stats.csv"), index = None)

#%%
#
#typ = "edge"
#
#if typ == "cell":
#    type_cell = yprof
#    dst = os.path.join(dst, "{}_guassian_fit_cell_profiles.pdf".format(typ))
#else:
#    type_cell = yprof_e
#    dst = os.path.join(dst, "{}_guassian_fit_edge_profiles.pdf".format(typ))
#
#
#pdf_pages = PdfPages(dst) #compiles into multiple pdfs
#
#
#for i in range(len(type_cell)):
#    try:
#        fig = plt.figure(figsize=(8,4))
#        ax = fig.add_axes([.4,.1,.5,.8])
#        
#        prof = type_cell[i][7:14]
#        x = ar(range(len(prof)))
#        y = prof
#        
#        # weighted arithmetic mean (corrected - check the section below)
#        mean = sum(x * y) / sum(y)
#        sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
#        
#        def Gauss(x, a, x0, sigma):
#            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
#        
#        popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma])
#        c = Gauss(x, *popt)
#        
#        textstr = "\n".join((
#                  "width (px): {:0.1f}".format(max(c)-min(c)),
#                  "left minima: {}".format(prof[0]),
#                  "right minima: {}".format(prof[len(prof)-1])))
#        
#        ax.plot(x, y, 'b+:', label='data')
#        ax.plot(x, c, 'r-', label='fit')
#        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#        # place a text box in upper left in axes coords
#        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
#        
#        ytick_spacing = 100; xtick_spacing = 1
#        ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
#        ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
#        ax.set_ylim([min(prof)-50, max(prof)+50])
#        ax.set_xlim([0, 6])
#        ax.set_xlabel('Distance (pixels)')
#        ax.set_ylabel('Intensity')
#        
#        pdf_pages.savefig(dpi = 300, bbox_inches = 'tight') 
#        plt.close()
#    except:
#        print(i, type_cell[i])
#
##write PDF document contains all points
#pdf_pages.close()
