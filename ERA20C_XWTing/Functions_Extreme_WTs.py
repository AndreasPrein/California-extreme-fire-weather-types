#!/usr/bin/env python
'''File name: Functions_Extreme-WTs.py
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 21.07.2019
    Date last modified: 21.07.2019

    ############################################################## 
    Purpos:

    Contains functions that are called for the Weather Typing
    of extreme precipiation events

'''

from dateutil import rrule
import datetime
import glob
from netCDF4 import Dataset
import sys, traceback
import dateutil.parser as dparser
import string
from pdb import set_trace as stop
import numpy as np
import numpy.ma as ma
import os
from mpl_toolkits import basemap
# import ESMF
import pickle
import subprocess
import pandas as pd
from scipy import stats
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import pylab as plt
import random
import scipy.ndimage as ndimage
import scipy
import shapefile
import matplotlib.path as mplPath
from matplotlib.patches import Polygon as Polygon2
# Cluster specific modules
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.vq import kmeans2,vq, whiten
from scipy.ndimage import gaussian_filter
import seaborn as sns
# import metpy.calc as mpcalc
import matplotlib.gridspec as gridspec
from collections import OrderedDict

import warnings
warnings.filterwarnings("ignore")

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height]) #,axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax



def Scatter_ED_PR(EuclDist,
                  ClosestWT,
                  PRall,
                  NrExtremes,
                  PlotLoc='./',
                  PlotName='Scatter.pdf'):
    import matplotlib.gridspec as gridspec
    from matplotlib import pyplot
    # plots a scatter diagram comparing eucledian distacnes from 
    # EWT centroids with daily precipitation volumes in the target region

    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize=(15,15), constrained_layout=True)
    widths = [3, 1]
    heights = [1, 3]
    gs1 = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=widths, height_ratios=heights)

    # plot scatter first
    ax = fig.add_subplot(gs1[1, 0])
    ax.scatter(EuclDist,PRall, color="k", s=2)
    # highlight the top extreme events
    ExtremePR=PRall[np.argsort(PRall)][::-1][:NrExtremes]
    ExtrPR_ED=EuclDist[np.argsort(PRall)][::-1][:NrExtremes]
    # lable the WTs of the most extreme days in color
    WTcolors=['#1f78b4','#33a02c','#e31a1c','#ff7f00','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#6a3d9a','#ffff99','#b15928','#1f78b4','#33a02c','#e31a1c','#ff7f00','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#6a3d9a','#ffff99','#b15928','#1f78b4','#33a02c','#e31a1c','#ff7f00','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#6a3d9a','#ffff99','#b15928','#1f78b4','#33a02c','#e31a1c','#ff7f00','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#6a3d9a','#ffff99','#b15928']
    XWTs=ClosestWT[np.argsort(PRall)][::-1][:NrExtremes]
    XWTunique=np.unique(XWTs)
    for xwt in range(len(XWTunique)):
        iAct=(XWTs == XWTunique[xwt])
        try:
            ax.scatter(ExtrPR_ED[iAct], ExtremePR[iAct], alpha=1, color=WTcolors[XWTunique[xwt]], label='XWT-'+str(XWTunique[xwt]+1), s=30)
        except:
            stop()
    # ax.scatter(ExtrPR_ED,ExtremePR,color='r', s=15)
    Q75_Extr=np.percentile(ExtrPR_ED, 75)
    plt.axvline(x=Q75_Extr, c='r', ls='--', lw=2.5)

    plt.xlabel('Eucledian Distance []')
    plt.ylabel('SCS count [events per day]')
    Xrange=ax.get_xlim()
    Yrange=ax.get_ylim()
    ax.legend(loc='upper right')

    # plot histogram for EDs
    ax = fig.add_subplot(gs1[0, 0])
    bins = np.linspace(Xrange[0], Xrange[1], 50)
    pyplot.hist(EuclDist, bins, alpha=0.5, color='k',label='all events', density=True)
    pyplot.hist(ExtrPR_ED, bins, alpha=0.5, color='r', label='extremes', density=True)
    pyplot.legend(loc='upper right')
    plt.ylabel('Probability []')
    plt.xlim([Xrange[0],Xrange[1]])

    # plot histogram for PR
    ax = fig.add_subplot(gs1[1, 1])
    bins = np.linspace(Yrange[0], Yrange[1], 50)
    plt.hist(np.array(PRall,dtype=np.float32), bins, alpha=0.5, color='k',label='all events', density=True, orientation="horizontal",histtype='step')
    pyplot.hist(np.array(ExtremePR, dtype=np.float32), bins, alpha=0.5, color='r', label='extremes', density=True, orientation="horizontal",histtype='step')
    plt.xlabel('Probability []')
    plt.ylim([Yrange[0],Yrange[1]])


    # Calculate the skill scores
    # from Functions_Extreme_WTs import ExtremeDays
    # rgiExtrEval=ExtremeDays(PRall,NrExtremes,7)
    rgiExtrEval=np.argsort(PRall)[-NrExtremes:]
    from Functions_Extreme_WTs import MRR, MRD, perkins_skill
    # Perkins Skill Score
    rPSS=perkins_skill(EuclDist,EuclDist[rgiExtrEval], 0.5)
    # Mean relative difference
    rMRD=MRD(EuclDist,PRall,rgiExtrEval)
    # Mean Rank Ratio
    rMRR=MRR(EuclDist,rgiExtrEval)
    # % of days excluded
    Excluded=(1-np.sum(EuclDist < np.nanpercentile(EuclDist[rgiExtrEval],75))/float(len(EuclDist)))*100.
    
    
    
#     # Mean relative difference
#     grMRD=MRD(EuclDist,testing_predictant,rgiExtrEval)
#     # Mean Rank Ratio
#     grMRR=MRR(EuclDist,rgiExtrEval)
#     # % of days excluded
#     grExluded=(1-np.sum(MinDistance < np.nanpercentile(EuclDist[rgiExtrEval],75))/float(len(EuclDist)))*100.
    # calculate the AUC
    from sklearn.metrics import roc_auc_score
    testy=(PRall >= np.sort(PRall)[-NrExtremes])
    probs=(EuclDist-np.min(EuclDist)); probs=np.abs((probs/probs.max())-1)
    try:
        auc = roc_auc_score(testy, probs)
    except:
        auc=np.nan

    # Calculate the Average precision-recall score
    from sklearn.metrics import average_precision_score
    from sklearn import svm, datasets
    try:
        average_precision = average_precision_score(testy, probs)
    except:
        average_precision = np.nan
    

    # Write skill-scores in top right corner of plot
    plt.text(0.70, 0.85, 'Perkins Skill Score: '+str("%.2f" % rPSS), fontsize=17, transform=plt.gcf().transFigure)
    plt.text(0.70, 0.81, 'Mean relative difference: '+str("%.2f" % rMRD), fontsize=17, transform=plt.gcf().transFigure)
    plt.text(0.70, 0.77, 'Mean Rank Ratio: '+str("%.2f" % rMRR), fontsize=17, transform=plt.gcf().transFigure)
    plt.text(0.70, 0.73, 'Excluded days: '+str("%.2f" % Excluded)+' %', fontsize=17, transform=plt.gcf().transFigure)
    plt.text(0.70, 0.69, 'AUC: '+str("%.2f" % auc)+' %', fontsize=17, transform=plt.gcf().transFigure)
    plt.text(0.70, 0.65, 'APR: '+str("%.2f" % average_precision)+' %', fontsize=17, transform=plt.gcf().transFigure)
    print('    plot: '+PlotLoc+PlotName)
    fig.savefig(PlotLoc+PlotName)



def ReadPRISM(rgiYears,        # array containing the years that should be read
              iNrOfExtremes,   # number of extreme events
              rgiSeasonWT,     # months that should be processed
              iMonths,         # array of months that should be read in
              grPRregion):   # shapefile that contains target region

    MinDistDD = 0
    ncid=Dataset('/glade/campaign/mmm/c3we/prein/observations/PRISM/data/PR/PRISM_daily_ppt_2014.nc', mode='r') # open the netcdf file
    rgrLatPR=np.squeeze(ncid.variables['lat'][:])
    rgrLonPR=np.squeeze(ncid.variables['lon'][:])
    ncid.close()
    rgrGridCells=[(rgrLonPR.ravel()[ii],rgrLatPR.ravel()[ii]) for ii in range(len(rgrLonPR.ravel()))]
    rgrSRact=np.array(grPRregion.contains_points(rgrGridCells)); rgrSRact=np.reshape(rgrSRact, (rgrLatPR.shape[0], rgrLatPR.shape[1]))
    rgiSrPR=np.array(np.where(rgrSRact == True))
    iLatMax=rgiSrPR[0,:].max()
    iLatMin=rgiSrPR[0,:].min()
    iLonMax=rgiSrPR[1,:].max()
    iLonMin=rgiSrPR[1,:].min()
    rgrPRdata=np.zeros((sum(rgiSeasonWT),iLatMax-iLatMin,iLonMax-iLonMin))
    jj=0
    for yy in range(len(rgiYears)):
        rgdTimeYY = pd.date_range(datetime.datetime(rgiYears[0]+yy, 1, 1,0), end=datetime.datetime(rgiYears[0]+yy, 12, 31,23), freq='d')
        rgiDD=np.where(((rgdTimeYY.year == rgiYears[0]+yy) & (np.isin(rgdTimeYY.month, iMonths))))[0]
        ncid=Dataset('/glade/campaign/mmm/c3we/prein/observations/PRISM/data/PR/PRISM_daily_ppt_'+str(rgiYears[0]+yy)+'.nc', mode='r') # open the netcdf file
        rgrPRdata[jj:jj+len(rgiDD),:,:]=np.squeeze(ncid.variables['PR'][rgiDD,iLatMin:iLatMax,iLonMin:iLonMax])
        ncid.close()
        jj=jj+len(rgiDD)
    rgrPRdata[rgrPRdata < 0]=np.nan
    rgiSRgridcells=rgrSRact[iLatMin:iLatMax,iLonMin:iLonMax]
    
    if iNrOfExtremes >= rgrPRdata.shape[0]:
        iNrOfExtremes = rgrPRdata.shape[0]

    rgrPRrecords=np.nanmean(rgrPRdata[:,rgiSRgridcells], axis=(1))
    SortedDates=np.argsort(rgrPRrecords)[:][::-1]
    rgiExtremePR=np.zeros((iNrOfExtremes)); rgiExtremePR[:]=np.nan
    ii=1
    jj=1
    rgiExtremePR[0]=SortedDates[0]
    while ii < iNrOfExtremes:
        if np.nanmin(np.abs(rgiExtremePR - SortedDates[jj])) < MinDistDD:
            jj=jj+1
        else:
            rgiExtremePR[ii]=SortedDates[jj]
            jj=jj+1
            ii=ii+1
    rgiExtremePR=rgiExtremePR.astype('int')

    return rgrPRrecords, rgiExtremePR, rgrPRdata[:,rgiSRgridcells]





# ###################################################
def ReadERAI(grWTregion,
             rgdTime,
             iMonths,         # list of months that should be considered
             rgsWTfolders,
             rgsWTvars):

    rgiYears=np.unique(rgdTime.year)

    ncid=Dataset('/glade/scratch/prein/ERA-Interim/PSL/fin_PSL-sfc_ERA-Interim_12-0_2014.nc', mode='r') # open the netcdf file
    rgrLatWT1D=np.squeeze(ncid.variables['lat'][:])
    rgrLonWT1D=np.squeeze(ncid.variables['lon'][:])
    ncid.close()
    rgrLonWT=np.asarray(([rgrLonWT1D,]*rgrLatWT1D.shape[0]))
    rgrLonWT[rgrLonWT > 180]=rgrLonWT[rgrLonWT > 180]-360
    rgrLatWT=np.asarray(([rgrLatWT1D,]*rgrLonWT1D.shape[0])).transpose()
    
    rgrGridCells=[(rgrLonWT.ravel()[ii],rgrLatWT.ravel()[ii]) for ii in range(len(rgrLonWT.ravel()))]
    rgrSRact=np.array(grWTregion.contains_points(rgrGridCells)); rgrSRact=np.reshape(rgrSRact, (rgrLatWT.shape[0], rgrLatWT.shape[1]))
    rgiSrWT=np.array(np.where(rgrSRact == True))
    iLatMax=rgiSrWT[0,:].max()
    iLatMin=rgiSrWT[0,:].min()
    iLonMax=rgiSrWT[1,:].max()
    iLonMin=rgiSrWT[1,:].min()
    
    DailyVars=np.zeros((len(rgdTime),iLatMax-iLatMin,iLonMax-iLonMin,len(rgsWTvars))); DailyVars[:]=np.nan
    for yy in range(len(rgiYears)):
        print('        Read ERA-I year: '+str(rgiYears[yy]))
        DaysYY = pd.date_range(datetime.datetime(rgiYears[yy], 1, 1,0), end=datetime.datetime(rgiYears[yy], 12, 31,23), freq='d')
        DD=((rgdTime.year == rgiYears[yy]) & np.isin(rgdTime.month, iMonths))
        DDactYYYY=np.isin(DaysYY.month, iMonths)
        # DDactYYYY=((DaysYY.month >= iStartMon) & (DaysYY.month <= iStopMon))
        for va in range(len(rgsWTvars)):
            ncid=Dataset(rgsWTfolders[va]+str(rgiYears[yy])+'.nc', mode='r')
            try:
                DailyVars[DD,:,:,va]=np.squeeze(np.squeeze(ncid.variables[rgsWTvars[va]])[:,iLatMin:iLatMax,iLonMin:iLonMax])[DDactYYYY,:]
            except:
                stop()
            ncid.close()


    return DailyVars, rgrLonWT[iLatMin:iLatMax,iLonMin:iLonMax], rgrLatWT[iLatMin:iLatMax,iLonMin:iLonMax]

# ###################################################
def ReadERA5(grWTregion,
             rgdTime,
             iMonths,         # list of months that should be considered
             rgsWTfolders,
             rgsWTvars):

    rgiYears=np.unique(rgdTime.year)

    ncid=Dataset('/glade/campaign/mmm/c3we/prein/ERA5/U850/U850_2006.nc', mode='r') # open the netcdf file
    rgrLatWT1D=np.squeeze(ncid.variables['latitude'][:])
    rgrLonWT1D=np.squeeze(ncid.variables['longitude'][:])
    ncid.close()
    rgrLonWT=np.asarray(([rgrLonWT1D,]*rgrLatWT1D.shape[0]))
    rgrLonWT[rgrLonWT > 180]=rgrLonWT[rgrLonWT > 180]-360
    rgrLatWT=np.asarray(([rgrLatWT1D,]*rgrLonWT1D.shape[0])).transpose()
    
    rgrGridCells=[(rgrLonWT.ravel()[ii],rgrLatWT.ravel()[ii]) for ii in range(len(rgrLonWT.ravel()))]
    rgrSRact=np.array(grWTregion.contains_points(rgrGridCells)); rgrSRact=np.reshape(rgrSRact, (rgrLatWT.shape[0], rgrLatWT.shape[1]))
    rgiSrWT=np.array(np.where(rgrSRact == True))
    iLatMax=rgiSrWT[0,:].max()
    iLatMin=rgiSrWT[0,:].min()
    iLonMax=rgiSrWT[1,:].max()
    iLonMin=rgiSrWT[1,:].min()
    
    DailyVars=np.zeros((len(rgdTime),iLatMax-iLatMin,iLonMax-iLonMin,len(rgsWTvars))); DailyVars[:]=np.nan
    for yy in range(len(rgiYears)):
        print('        Read ERA-5 year: '+str(rgiYears[yy]))
        DaysYY = pd.date_range(datetime.datetime(rgiYears[yy], 1, 1,0), end=datetime.datetime(rgiYears[yy], 12, 31,23), freq='d')
        DD=((rgdTime.year == rgiYears[yy]) & np.isin(rgdTime.month, iMonths))
        DDactYYYY=np.isin(DaysYY.month, iMonths)
        # DDactYYYY=((DaysYY.month >= iStartMon) & (DaysYY.month <= iStopMon))
        for va in range(len(rgsWTvars)):
            ncid=Dataset(rgsWTfolders[va]+rgsWTvars[va]+'_'+str(rgiYears[yy])+'.nc', mode='r')
            try:
                DailyVars[DD,:,:,va]=np.squeeze(np.squeeze(ncid.variables[rgsWTvars[va]])[:,iLatMin:iLatMax,iLonMin:iLonMax])[DDactYYYY,:]
            except:
                TMP=np.squeeze(np.squeeze(ncid.variables[rgsWTvars[va]])[:,iLatMin:iLatMax,iLonMin:iLonMax])[DDactYYYY,:]
                DailyVars[DD,:,:,va]=TMP[:DailyVars[DD,:,:,va].shape[0],:]
            ncid.close()
            if (rgsWTvars[va] == 'CAPE') | (rgsWTvars[va] == 'CIN'):
                CACI = DailyVars[DD,:,:,va]
                CACI[np.isnan(CACI)] = 0
                DailyVars[DD,:,:,va] = CACI

    return DailyVars, rgrLonWT[iLatMin:iLatMax,iLonMin:iLonMax], rgrLatWT[iLatMin:iLatMax,iLonMin:iLonMax]


# ###################################################
def PreprocessWTdata(DailyVarsInput,                      # WT data [time,lat,lon,var]
                     RelAnnom=1,                     # calculate relative anomalies [1-yes; 0-no]
                     SmoothSigma=0,                  # Smoothing stddev (Gaussian smoothing)
                     RemoveAnnualCycl=1,             # remove annual cycle [1-yes; 0-no]
                     NormalizeData='D',              # normalize variables | options are  - 'C' - climatology
                                                     # - 'D' - daily (default)
                                                     # - 'N' - none
                     ReferencePer=None):             # period for normalizing the data
    

    DailyVars = np.copy(DailyVarsInput)
    # Calculate relative anomaly
    if RelAnnom == 1:
        # we have to work with absolute values for this since we risk to divide by zero values in the climatology
        DailyVars=np.abs(DailyVars)
        if ReferencePer is None:
            DailyVars=(DailyVars-np.mean(DailyVars, axis=0)[None,:])/np.mean(DailyVars, axis=0)[None,:]
        else:
            DailyVars=(DailyVars-np.mean(DailyVars[ReferencePer], axis=0)[None,:])/np.mean(DailyVars[ReferencePer], axis=0)[None,:]

    if len(DailyVars.shape) == 3:
        DailyVars = DailyVars[:,:,:,None]
    # Spatially smooth the data
    DailyVars=gaussian_filter(DailyVars[:,:,:,:], sigma=(0,SmoothSigma,SmoothSigma,0))

    # Remove the annual cycle
    if RemoveAnnualCycl == 1:
        SpatialMeanData=pd.DataFrame(np.nanmean(DailyVars, axis=(1,2)))
        Averaged=np.roll(np.array(SpatialMeanData.rolling(window=21).mean()), -10, axis=0)
        Averaged[:10,:]=Averaged[11,:][None,:]; Averaged[-10:,:]=Averaged[-11,:][None,:]
        DailyVars=DailyVars-Averaged[:,None,None,:]

    # Normalize the data
    if NormalizeData == 'D':
        DailyVars=(DailyVars-np.mean(DailyVars, axis=(1,2))[:,None,None,:])/np.std(DailyVars, axis=(1,2))[:,None,None,:]
    elif NormalizeData == 'C':
        if ReferencePer is None:
            DailyVars=(DailyVars-np.mean(DailyVars, axis=(0,1,2))[None,None,None,:])/np.std(DailyVars, axis=(0,1,2))[None,None,None,:]
        else:
            DailyVars=(DailyVars-np.mean(DailyVars[ReferencePer], axis=(0,1,2))[None,None,None,:])/np.std(DailyVars[ReferencePer], axis=(0,1,2))[None,None,None,:]
        DailyVars[np.isnan(DailyVars)]=0

    return DailyVars




# ===================================================================
def GetExtremeDays(DailyVars,
                   rgdTime,
                   rgiExtremeDays):
    # Grab the extreme days from the full data
    rgrWTdata=np.zeros((len(rgiExtremeDays),DailyVars.shape[1],DailyVars.shape[2],DailyVars.shape[3]))
    for dd in range(len(rgiExtremeDays)):
        rgdTimeYY = pd.date_range(datetime.datetime(rgiExtremeDays[dd].year, 1, 1,0), end=datetime.datetime(rgiExtremeDays[dd].year, 12, 31,23), freq='d')
        rgiDD=np.where(((rgdTime.year == rgiExtremeDays[dd].year) & (rgdTime.month ==rgiExtremeDays[dd].month ) & (rgdTime.day == rgiExtremeDays[dd].day)))[0]
        rgrWTdata[dd,:,:,:]=DailyVars[rgiDD[0],:]

    return rgrWTdata





# ===================================================================
def ClusterAnalysis(rgrWTdata,
                    sPlotDir,
                    iNrOfExtremes,
                    YYYY_stamp,
                    Plot=0,        # 0 - no plots; 1 - plots will be saved to sPlotDir
                    ClusterMeth='HandK',   # current options are ['HandK','hdbscan']
                    ClusterBreakup=0):

    if ClusterMeth == 'hdbscan':
    # #--------------------------------------------------
        # # HDBSCAN -- https://hdbscan.readthedocs.io/en/latest/basic_hdbscan.html
        import hdbscan
        rgrDataCluster=np.reshape(rgrWTdata, (rgrWTdata.shape[0], rgrWTdata.shape[1]*rgrWTdata.shape[2]*rgrWTdata.shape[3]))
        
        # clusterer = hdbscan.HDBSCAN(min_cluster_size=2).fit(rgrDataCluster) 
        try:
            Epsilon = int(np.std(np.sum(np.abs(rgrDataCluster), axis=1)))
        except:
            rgrClustersFin = None
            return rgrClustersFin
        MinClusterSize = 2
#         clusterer = hdbscan.HDBSCAN(min_cluster_size=MinClusterSize, min_samples=1, cluster_selection_epsilon=Epsilon).fit(rgrDataCluster)
#         clusterer = hdbscan.HDBSCAN().fit(rgrDataCluster)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=MinClusterSize).fit(rgrDataCluster)
        Clusters=clusterer.labels_

        # Recluster the outlier cases seperately
        Outlier = (Clusters == -1)
        if np.sum(Outlier) == len(Clusters):
            Clusters[:]=0
        
        if np.sum(Outlier) > 2:
            ClustersO = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1, cluster_selection_epsilon=1).fit(rgrDataCluster[Outlier,:]).labels_
            ClustersT = ClustersO + Clusters.max()+1
            ClustersT[ClustersO == -1] = -1
            Clusters[Outlier] = ClustersT
        
        if ClusterBreakup == 1:
            # split a cluster if it contains unproportionally many cases
            ClusterSizes = np.array([np.sum(Clusters == ii) for ii in range(np.max(Clusters))])
            if len(ClusterSizes) == 0:
                ClusterSizes = [len(Clusters)]
            try:
                MAX = max(ClusterSizes)
            except:
                stop()
            if (max(ClusterSizes)/len(Clusters) > 0.5) & (len(Clusters) > 10) == True:
                iClusterL = (Clusters == np.argmax(ClusterSizes))
                ClusterL = hdbscan.HDBSCAN(min_cluster_size=3, min_samples=1, cluster_selection_epsilon=1).fit(rgrDataCluster[iClusterL,:]).labels_
                # Recluster the outlier cases seperately
                Outlier = (ClusterL == -1)
            
                if np.sum(Outlier) > 2:
                    ClustersO = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1, cluster_selection_epsilon=Epsilon).fit(rgrDataCluster[iClusterL,:][Outlier,:]).labels_
                    ClustersO = ClustersO + ClusterL.max()+1
                    ClusterL[Outlier] = ClustersO
    
                ClusterL = ClusterL + Clusters.max()+1
                Clusters[iClusterL] = ClusterL
                Clusters[Clusters >= np.argmax(ClusterSizes)] = Clusters[Clusters >= np.argmax(ClusterSizes)]-1

        UniqueClus=np.where(Clusters == -1)[0]
        if len(UniqueClus) > 0:
            MaxClust=Clusters.max()
            for uc in range(len(UniqueClus)):
                Clusters[UniqueClus[uc]]=MaxClust+1+uc
        # check if cluster count starts with zero
        if np.min(Clusters) != 0:
            Clusters=Clusters-np.min(Clusters)
    
        # calculate centroids
        rgrWTcentroids=np.zeros((Clusters.max()+1, rgrWTdata.shape[1], rgrWTdata.shape[2], rgrWTdata.shape[3])); rgrWTcentroids[:]=np.nan
        for cc in range(Clusters.max()+1):
            rgiClAct=(Clusters == (cc))
            rgrWTcentroids[cc,:]=np.mean(rgrWTdata[rgiClAct,:,:,:], axis=0)
        rgrWTcentroids=np.reshape(rgrWTcentroids, (rgrWTcentroids.shape[0], rgrWTcentroids.shape[1]*rgrWTcentroids.shape[2]*rgrWTcentroids.shape[3]))
        rgrClustersFin=(rgrWTcentroids, Clusters, rgrWTdata)

    if ClusterMeth == 'HandK':
        # --------------------------------------------------
        # # PERFORM HIRACHICAL AND K-MEANS CLUSTER ANALYSIS
        # see excample: https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
        rgrDataCluster=np.reshape(rgrWTdata, (rgrWTdata.shape[0],rgrWTdata.shape[1]*rgrWTdata.shape[2]*rgrWTdata.shape[3]))
        rgrCluster = linkage(rgrDataCluster, 'ward')
        cc, coph_dists = cophenet(rgrCluster, pdist(rgrDataCluster))
        
        # last = rgrCluster[-25:, 2]
        last = rgrCluster[-6:, 2]
        last_rev = last[::-1]
        idxs = np.arange(1, len(last) + 1)
        
        acceleration = np.diff(last, 2)  # 2nd derivative of the distances
        acceleration_rev = acceleration[::-1]
    
    
        if Plot == 1:
            # PLOT DENDROGRAM
            fig = plt.figure(figsize=(10, 6))
            plt.title('Hierarchical Clustering Dendrogram')
            plt.xlabel('sample index')
            plt.ylabel('distance')
            from Functions_Extreme_WTs import fancy_dendrogram
            fancy_dendrogram(rgrCluster,
                             truncate_mode='lastp',  # show only the last p merged clusters
                             p=12,  # show only the last p merged clusters
                             leaf_rotation=90.,  # rotates the x axis labels
                             leaf_font_size=12.,  # font size for the x axis labels
                             show_contracted=True,  # to get a distribution impression in truncated branches
                             annotate_above=10,
                             max_d=36,
                         )
            sPlotFile=sPlotDir
            sPlotName= 'BottomUp-Hirarch-Cluster_Dendrogram_'+str("%03d" % iNrOfExtremes)+'_'+YYYY_stamp+'.pdf'
            if os.path.isdir(sPlotFile) != 1:
                subprocess.call(["mkdir","-p",sPlotFile])
            print('        Plot map to: '+sPlotFile+sPlotName)
            fig.savefig(sPlotFile+sPlotName)
        
            # PLOT DISTANCE ACCELERATION
            fig = plt.figure(figsize=(7, 6))
            plt.title('Hierarchical Clustering Dendrogram')
            plt.xlabel('Nr. of weather types')
            plt.ylabel('distance')
            plt.plot(idxs, last_rev, c='k')
    
            plt.plot(idxs[:-2] + 1, acceleration_rev, c='r')
            sPlotFile=sPlotDir
            sPlotName= 'Accelaration-Of-Distance-Growth_'+str("%03d" % iNrOfExtremes)+'_'+YYYY_stamp+'.pdf'
            if os.path.isdir(sPlotFile) != 1:
                subprocess.call(["mkdir","-p",sPlotFile])
            print('        Plot map to: '+sPlotFile+sPlotName)
            fig.savefig(sPlotFile+sPlotName)
        iClusters1=0
        try:
            while acceleration_rev[iClusters1] > 0: iClusters1=iClusters1+1
            iClusters2=np.where(max(acceleration_rev) == acceleration_rev)[0][0]+1
        except:
            iClusters2=iClusters1
        iClusters=np.max([iClusters1,iClusters2])+1
        iClusters=3 #np.max([iClusters1,iClusters2])+1
        rThreshold=last[-iClusters]
    
        rgiClusterMembers=fcluster(rgrCluster,rThreshold, criterion='distance')
        rgrClusters=np.zeros((rgiClusterMembers.max(),rgrDataCluster.shape[1])); rgrClusters[:]=np.nan
        for cc in range(rgiClusterMembers.max()):
            rgiClAct=(rgiClusterMembers == (cc+1))
            rgrClusters[cc,:]=np.mean(rgrDataCluster[rgiClAct,:], axis=0)

        # use this as initial seed for the k-means clustering
        rgrClustersFin=kmeans2(rgrDataCluster,rgrClusters)
    return rgrClustersFin





# ===================================================================
def EucledianDistance(DailyVars,
                      rgrClustersFin):
    SHAPE=DailyVars.shape
    Data_flatten=np.reshape(DailyVars, (SHAPE[0],SHAPE[1]*SHAPE[2]*SHAPE[3]))
    EucledianDist=np.zeros((SHAPE[0],rgrClustersFin[0].shape[0])); EucledianDist[:]=np.nan
    Correlation=np.copy(EucledianDist)
    for dd in range(SHAPE[0]):
        EucledianDist[dd,:] = np.array([np.linalg.norm(rgrClustersFin[0][wt,:]-Data_flatten[dd,:]) for wt in range(rgrClustersFin[0].shape[0])])
        Correlation[dd,:] = np.array([np.corrcoef(rgrClustersFin[0][wt,:],Data_flatten[dd,:])[0][1] for wt in range(rgrClustersFin[0].shape[0])])
    
    return EucledianDist, Correlation




# Perkins Skill Score
# overlap between two PDFs | zero is best
def perkins_skill(data1, data2, Binsize):
    Min=np.nanmin([np.nanmin(data1),np.nanmin(data2)])
    Max=np.nanmax([np.nanmax(data1),np.nanmax(data2)])

    hist, bin_edges = np.histogram(data1[~np.isnan(data1)],bins=np.arange(Min,Max,Binsize),density=True)
    pdf1 = hist*np.diff(bin_edges)
    try:
        histEx, bin_edgesEx = np.histogram(data2,bins=np.arange(Min,Max,Binsize),density=True)
    except:
        stop()
    pdf2 = histEx*np.diff(bin_edgesEx)
    mins = np.minimum(pdf1,pdf2)
    ss = np.nansum(mins)
    return ss

# Mean relative difference
# Mean relative difference between extreme PR cases and cases with low
# Eucledian Distances | zero is best
def MRD(Distance, Precipitation, iExtremes):
    Extreme75P=np.nanpercentile(Distance[iExtremes],75)
    MeanExtreme=np.nanmean(np.array(Precipitation[iExtremes][Distance[iExtremes] <= Extreme75P],dtype=np.float32))
    MeanAll75P=np.nanmean(np.array(Precipitation[Distance <= Extreme75P],dtype=np.float32))
    RelDiff=((MeanAll75P-MeanExtreme)/MeanExtreme)*-1
    return RelDiff

# Mean Rank Ratio
# Difference between ranks of extreme cases according to their Eucledian Distances
# and average rank in dataset | zero is best; one is no scill; 2 is perfect negative skill
def MRR(Distance,iExtremes):
    RankedDistances=np.argsort(Distance)
    ExtremeRanks=np.mean(np.array([np.where(RankedDistances == iExtremes[ii])[0][0] for ii in range(len(iExtremes))]))-(len(iExtremes)/2.-0.5)
    WorstRanks=len(Distance)-1
    MeanRankRatio=(ExtremeRanks/(WorstRanks-len(iExtremes)+1))*2
    return MeanRankRatio

# Get extreme PR days with minimum days appart
def ExtremeDays(Record,ExtremeNr, DistanceDD):
    SortedDates=np.argsort(Record)[:][::-1]
    rgiExtremePR=np.zeros((ExtremeNr)); rgiExtremePR[:]=np.nan
    ii=1
    jj=1
    rgiExtremePR[0]=SortedDates[0]
    while ii < ExtremeNr:
        if np.nanmin(np.abs(rgiExtremePR - SortedDates[jj])) < DistanceDD:
            jj=jj+1
        else:
            rgiExtremePR[ii]=SortedDates[jj]
            jj=jj+1
            ii=ii+1
    return rgiExtremePR.astype('int')




# ===================================================================
def read_shapefile(sf):
    """
    Read a shapefile into a Pandas dataframe with a 'coords' 
    column holding the geometry information. This uses the pyshp
    package
    """
    fields = [x[0] for x in sf.fields][1:]
    records = sf.records()
    shps = [s.points for s in sf.shapes()]
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    return df





# ===================================================================
def ReadCESMday(DaySel,
                Exp,
                iWest,
                iEast,
                iSouth,
                iNort,
                rgrTimeCESMFull,
                VARS=None,
                AddCells=0):

    """
    Read in a single day within a region from one
    CESM large ensemble simulation
    All variables nescessary for a synopic mapplot are read in
    """

    if VARS == None:
        rgsWTvars=['Z500','U850','V850','TMQ',]
        VarsFullName=  ['Z500','U850','V850','PW']
        rgsWTfolders=['/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/Z500/',\
                      '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/U850/',\
                      '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/V850/',\
                      '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/TMQ/']
    else:
        rgsWTvars=VARS[0]
        VarsFullName=VARS[1]
        rgsWTfolders=VARS[2]
    s20Cname='b.e11.B20TRC5CNBDRD.f09_g16.'
    s21Cname='b.e11.BRCP85C5CNBDRD.f09_g16.'

    # start reading in the CESM data
    iRegionPlus=AddCells # grid cell added around shape rectangle
    ncid=Dataset('/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/PSL/b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h1.PSL.18500101-20051231.nc', mode='r')
    rgrLonWT1D=np.squeeze(ncid.variables['lon'][:])
    rgrLatWT1D=np.squeeze(ncid.variables['lat'][:])
    ncid.close()
    rgrLonS=rgrLonWT1D[iWest-iRegionPlus:iEast+iRegionPlus]
    rgrLatS=rgrLatWT1D[iSouth-iRegionPlus:iNort+iRegionPlus]

    # Read the variables
    DataAll=np.zeros((len(rgrLatS),len(rgrLonS),len(rgsWTvars))); DataAll[:]=np.nan
    for va in range(len(rgsWTvars)):
        if DaySel.year < 2006:
            if Exp  == '001':
                rgrTimeCESM=pd.date_range(datetime.date(1850, 1, 1), end=datetime.date(2005, 12, 31), freq='d')
            else:
                rgrTimeCESM=pd.date_range(datetime.date(1920, 1, 1), end=datetime.date(2005, 12, 31), freq='d')
            Cfiles=glob.glob(rgsWTfolders[va]+'/'+s20Cname+Exp+'*'+rgsWTvars[va]+'*')[0]
        if DaySel.year >= 2006:
            if int(Exp)  >= 34:
                rgrTimeCESM=pd.date_range(datetime.date(2006, 1, 1), end=datetime.date(2100, 12, 31), freq='d')
                Cfiles=glob.glob(rgsWTfolders[va]+'/'+s21Cname+Exp+'*'+rgsWTvars[va]+'*')[0]
            elif DaySel.year <= 2080:
                rgrTimeCESM=pd.date_range(datetime.date(2006, 1, 1), end=datetime.date(2080, 12, 31), freq='d')
                try:
                    Cfiles=np.sort(glob.glob(rgsWTfolders[va]+'/'+s21Cname+Exp+'*'+rgsWTvars[va]+'*'))[0]
                except:
                    stop()
            elif DaySel.year >= 2081:
                rgrTimeCESM=pd.date_range(datetime.date(2081, 1, 1), end=datetime.date(2100, 12, 31), freq='d')
                Cfiles=np.sort(glob.glob(rgsWTfolders[va]+'/'+s21Cname+Exp+'*'+rgsWTvars[va]+'*'))[1]
        rgiNonLeap=np.where((rgrTimeCESM.month != 2) | (rgrTimeCESM.day != 29))[0]
        rgrTimeCESM=rgrTimeCESM[rgiNonLeap]
        iDDselect=np.where(rgrTimeCESM == DaySel)[0][0]

        
        try:
            ncid=Dataset(Cfiles, mode='r')
            DataAll[:,:,va]=np.squeeze(ncid.variables[rgsWTvars[va]][iDDselect,iSouth-iRegionPlus:iNort+iRegionPlus,iWest-iRegionPlus:iEast+iRegionPlus])
            ncid.close()
        except:
            stop()
        

    return DataAll, rgrLonS, rgrLatS


# ===================================================================
def DetrentData(DATA,
                TIME,
                YYYY_WINDOW):
    # remove YYYY_WINDOW anomaly from each year in the time series
    # This removes the thermodynamic effects while maintaining the dynamics effects
    rgrDataDetrended=np.copy(DATA)
    iYearsFull=np.unique(TIME.year)
    Y_half=int(YYYY_WINDOW)
    for yy in range(len(iYearsFull)):
        if yy < Y_half:
            yyStart=iYearsFull[0]
        else:
            yyStart=iYearsFull[yy-Y_half]
        if yy > len(iYearsFull)-(Y_half+1):
            yyStop=iYearsFull[len(iYearsFull)-1]
        else:
            yyStop=iYearsFull[yy+Y_half]
        iTimePeriodAct=((TIME.year >= yyStart) & (TIME.year <= yyStop))
        rgrDataDetrended[iTimePeriodAct,:,:,:]=DATA[iTimePeriodAct,:,:,:]-np.mean(DATA[iTimePeriodAct,:,:,:], axis=(0,1,2))[None,None,None,:]

        rgrDataDetrended=rgrDataDetrended+np.mean(DATA, axis=0)[None,:,:,:]

        return rgrDataDetrended


# ===================================================================
def SearchOptimum_XWT(PlotFile,
                      sPlotDir,
                      SkillScores_All,
                      GlobalMinimum1,
                      GlobalMinimum2,
                      Optimum,
                      VariableIndices,
                      Dimensions,
                      Metrics,
                      VarsFullName,
                      ss,
                      rgrNrOfExtremes,
                      WT_Domains,
                      Annual_Cycle,
                      SpatialSmoothing):
    # provide visual guidance for how large the impacts of changing
    # setup variables are to find an optimal configuration

    fig = plt.figure(figsize=(18, 13))
    plt.rcParams.update({'font.size': 14})
    rgsLableABC=list(string.ascii_lowercase)+list(string.ascii_uppercase)

    YY=[0,0,0,1,1]
    XX=[0,1,2,0,1]
    Shape=SkillScores_All.shape
    SkillScoreColors=['#1f78b4','#e31a1c','#ff7f00']

    # ---------------------------
    # Plot showing the time that extremes occured
    gs1 = gridspec.GridSpec(ncols=3, nrows=2, figure=fig)
    gs1.update(left=0.08, right=0.98,
               bottom=0.10, top=0.96,
               wspace=0.25, hspace=0.25)

    for se in range(len(Dimensions)-2):
        ax = fig.add_subplot(gs1[YY[se],XX[se]])
        if Dimensions[se] == 'Variables':
            SS1_Data=SkillScores_All[:,GlobalMinimum1[1][0],GlobalMinimum1[2][0],GlobalMinimum1[3][0],GlobalMinimum1[4][0],0,:]
            Xaxis1=[VarsFullName[int(VariableIndices[va,GlobalMinimum1[1][0],GlobalMinimum1[2][0],GlobalMinimum1[3][0],GlobalMinimum1[4][0],ss])] for va in range(Shape[0])]
            SS2_Data=SkillScores_All[:,GlobalMinimum2[1][0],GlobalMinimum2[2][0],GlobalMinimum2[3][0],GlobalMinimum2[4][0],1,:]
            Xaxis2=[VarsFullName[int(VariableIndices[va,GlobalMinimum2[1][0],GlobalMinimum2[2][0],GlobalMinimum2[3][0],GlobalMinimum2[4][0],ss])] for va in range(Shape[0])]
            Xlabel='variable'
            Xaxis=[Xaxis1[ii]+'/'+Xaxis2[ii] for ii in range(len(Xaxis1))]
            XOptimum=Optimum[0][0]
            YOptimum=SkillScores_All[Optimum[0][0],Optimum[1][0],Optimum[2][0],Optimum[3][0],Optimum[4][0],ss,:]
        if Dimensions[se] == 'Extreme Nr.':
            SS1_Data=SkillScores_All[GlobalMinimum1[0][0],:,GlobalMinimum1[2][0],GlobalMinimum1[3][0],GlobalMinimum1[4][0],0,:]
            SS2_Data=SkillScores_All[GlobalMinimum2[0][0],:,GlobalMinimum2[2][0],GlobalMinimum2[3][0],GlobalMinimum2[4][0],1,:]
            Xlabel='Nr. of extreme days'
            XOptimum=Optimum[1][0]
            YOptimum=SkillScores_All[Optimum[0][0],Optimum[1][0],Optimum[2][0],Optimum[3][0],Optimum[4][0],ss,:]
            Xaxis=rgrNrOfExtremes
        if Dimensions[se] == 'Domain Size':
            SS1_Data=SkillScores_All[GlobalMinimum1[0][0],GlobalMinimum1[1][0],:,GlobalMinimum1[3][0],GlobalMinimum1[4][0],0,:]
            SS2_Data=SkillScores_All[GlobalMinimum2[0][0],GlobalMinimum2[1][0],:,GlobalMinimum2[3][0],GlobalMinimum2[4][0],1,:]
            Xlabel='Domain size'
            XOptimum=Optimum[2][0]
            YOptimum=SkillScores_All[Optimum[0][0],Optimum[1][0],Optimum[2][0],Optimum[3][0],Optimum[4][0],ss,:]
            Xaxis=WT_Domains
        if Dimensions[se] == 'Annual Cycle':
            SS1_Data=SkillScores_All[GlobalMinimum1[0][0],GlobalMinimum1[1][0],GlobalMinimum1[2][0],:,GlobalMinimum1[4][0],0,:]
            SS2_Data=SkillScores_All[GlobalMinimum2[0][0],GlobalMinimum2[1][0],GlobalMinimum2[2][0],:,GlobalMinimum2[4][0],1,:]
            Xlabel='Annual cycle removed'
            XOptimum=Optimum[3][0]
            YOptimum=SkillScores_All[Optimum[0][0],Optimum[1][0],Optimum[2][0],Optimum[3][0],Optimum[4][0],ss,:]
            Xaxis=Annual_Cycle
        if Dimensions[se] == 'Smoothing':
            SS1_Data=SkillScores_All[GlobalMinimum1[0][0],GlobalMinimum1[1][0],GlobalMinimum1[2][0],GlobalMinimum1[3][0],:,0,:]
            SS2_Data=SkillScores_All[GlobalMinimum2[0][0],GlobalMinimum2[1][0],GlobalMinimum2[2][0],GlobalMinimum2[3][0],:,1,:]
            Xlabel='Spatial smoothing'
            XOptimum=Optimum[4][0]
            YOptimum=SkillScores_All[Optimum[0][0],Optimum[1][0],Optimum[2][0],Optimum[3][0],Optimum[4][0],ss,:]
            Xaxis=SpatialSmoothing

        # Start plotting the scill scores dependency on setup
        for sc in range(len(Metrics)):
            plt.plot(range(len(Xaxis)), SS1_Data[:,sc], c=SkillScoreColors[sc], label=Metrics[sc], ls='-')
            plt.plot(range(len(Xaxis)), SS2_Data[:,sc], c=SkillScoreColors[sc], label=Metrics[sc], ls='--')
            plt.plot(XOptimum, YOptimum[sc], marker='o', c=SkillScoreColors[sc], label='Optimum', markersize=10)
        plt.plot(range(len(Xaxis)), np.mean(SS1_Data[:,:], axis=1), c='k', label='Mean', ls='-',lw=3)
        plt.plot(range(len(Xaxis)), np.mean(SS2_Data[:,:], axis=1), c='k', label='Mean', ls='--',lw=3)
        plt.plot(XOptimum, np.mean(YOptimum), marker='o', c='k', label='Optimum', markersize=10)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel('Skill Score []')
        ax.set_xlabel(Xlabel)
        plt.xticks(np.arange(0, len(Xaxis), 1.0))
        ax.set_xticklabels(Xaxis, rotation=45)

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

    sPlotFile=sPlotDir
    # sPlotName= 'BottomUp-'+str(rgrClustersFin[1].max()+1)+'WT_precipitation.pdf'
    sPlotName= PlotFile
    if os.path.isdir(sPlotFile) != 1:
        subprocess.call(["mkdir","-p",sPlotFile])
    print('        Plot map to: '+sPlotFile+sPlotName)
    fig.savefig(sPlotFile+sPlotName)








# ===================================================================
def XWT(training_predictors, # predictor variables that are used to train the model
    testing_predictors,      # predictor variables that are used to evaluate the model
    training_predictant,     # predictent variable that is uesed to train the model
    testing_predictant,      # predictent variable that is uesed to evaluate the model
    training_time,           # daily time vector for the training dataset
    testing_time,            # daily time vector for the testing datast
    extreme_nr,              # Nr. of extreme events considered
    smoothing_radius,        # smoothing radius applied to predictor fields
    ClusterMeth='HandK',     # current options are ['HandK','hdbscan']
    ClusterBreakup=0,        # break up clusters that have more than 50% of events and if extreme Nr. > 10
    RelAnnom=1,              # 1 - calculate daily relative anomalies 
    NormalizeData=1,         # normalize daily variable fields - 1=yes
    MinDistDD=7,             # minimum nr of days between XWT events
    RemoveAnnualCycl=1):     # remove annual cycle in varaiables with 21 day moving average filter
    
    if len(training_predictant.shape) > 1:
        # data is given for grid cells
        training_predictant_GC = np.copy(training_predictant)
        testing_predictant_GC = np.copy(testing_predictant)
        training_predictant = np.nanmean(training_predictant_GC, axis=1)
        testing_predictant = np.nanmean(testing_predictant_GC, axis=1)

    #  OPTIONAL INPUTS
    sPlotDir=''
    YYYY_stamp=str(training_time[0].year)+'-'+str(training_time[0].year)

    from Functions_Extreme_WTs import ExtremeDays
    if extreme_nr == len(training_predictant):
        # all day weather typing
        rgiExtrTrain=np.array(range(len(training_predictant)))
        rgiExtrEval = np.array(range(len(testing_predictant)))
    else:
        # XWTing
        rgiExtrTrain=ExtremeDays(training_predictant,extreme_nr,MinDistDD)
        rgiExtrEval=np.argsort(testing_predictant)[-extreme_nr:]
    ExtrTrainDays=training_time[rgiExtrTrain]
    # rgiExtrEval=ExtremeDays(testing_predictant,extreme_nr,MinDistDD)

    from Functions_Extreme_WTs import PreprocessWTdata
    training_predictors=PreprocessWTdata(training_predictors,             # WT data [time,lat,lon,var]
                               RelAnnom=RelAnnom,                         # calculate relative anomalies [1-yes; 0-no]
                               SmoothSigma=smoothing_radius,              # Smoothing stddev (Gaussian smoothing)
                               RemoveAnnualCycl=RemoveAnnualCycl,         # remove annual cycle [1-yes; 0-no]
                               NormalizeData=NormalizeData)               # normalize data [1-yes; 0-no]

    from Functions_Extreme_WTs import GetExtremeDays
    try:
        rgrWTdata=GetExtremeDays(training_predictors,training_time,ExtrTrainDays)
    except:
        stop()
    
    # ################################################
    # ####  Run Hirarchical clustering
    from Functions_Extreme_WTs import ClusterAnalysis
    rgrClustersFin=ClusterAnalysis(rgrWTdata,
                                   sPlotDir,
                                   extreme_nr,
                                   YYYY_stamp,
                                   Plot=0,
                                   ClusterMeth=ClusterMeth,
                                   ClusterBreakup=ClusterBreakup)
    if rgrClustersFin == None:
        return None

    # ################################################
    # ####  Prepare evaluation data
    DailyVarsEvalNorm=PreprocessWTdata(testing_predictors,                  # WT data [time,lat,lon,var]
                                       RelAnnom=RelAnnom,                     # calculate relative anomalies [1-yes; 0-no]
                                       SmoothSigma=smoothing_radius,   # Smoothing stddev (Gaussian smoothing)
                                       RemoveAnnualCycl=RemoveAnnualCycl,             # remove annual cycle [1-yes; 0-no]
                                       NormalizeData=NormalizeData)                # normalize data [1-yes; 0-no]

    # ################################################
    # ######       EUCLEDIAN DISTANCES
    from Functions_Extreme_WTs import EucledianDistance
    if len(rgrClustersFin) == 3:
        # search clusters on centroids and individual XWT days
        rgrClustersFin = list(rgrClustersFin)
        rgrClustersFin[0] = np.append(rgrClustersFin[0],np.reshape(rgrClustersFin[2], (rgrClustersFin[2].shape[0],rgrClustersFin[2].shape[1]*rgrClustersFin[2].shape[2]*rgrClustersFin[2].shape[3])), axis=0)
        rgrClustersFin[1] = np.append(np.array(range(np.max(rgrClustersFin[1]+1))),rgrClustersFin[1])
        
       
    EucledianDist, Correlation =EucledianDistance(DailyVarsEvalNorm,
                                                  rgrClustersFin)

    from Functions_Extreme_WTs import Scatter_ED_PR
    MinDistance=np.min(EucledianDist, axis=1)
    ArgMIN = np.argmin(EucledianDist, axis=1)
    ClosestWT = [rgrClustersFin[1][ArgMIN[ii]] for ii in range(EucledianDist.shape[0])]
#     ClosestWT=np.argmin(EucledianDist, axis=1)
    MaxCorr=np.max(Correlation, axis=1)
    # Scatter_ED_PR(MinDistance,
    #               ClosestWT,
    #               Peval,
    #               rgrNrOfExtremes,
    #               PlotLoc=sPlotDir,
    #               PlotName='Scatter_'+sRegion+'_NrExt-'+str(rgrNrOfExtremes)+'_Smooth-'+str(SpatialSmoothing)+'_AnnCy-'+Annual_Cycle+'_'+VarsJoint+'_'+sMonths+'_'+Samples[ss]+'.pdf')

    # Calculate the skill scores
    from Functions_Extreme_WTs import MRR, MRD, perkins_skill
    # Perkins Skill Score
    try:
        grPSS=perkins_skill(MinDistance,MinDistance[rgiExtrEval], 0.5)
    except:
        grPSS=None
        
    # Mean relative difference
    grMRD=MRD(MinDistance,testing_predictant,rgiExtrEval)
    # Mean Rank Ratio
    grMRR=MRR(MinDistance,rgiExtrEval)
    # % of days excluded
    grExluded=(1-np.sum(MinDistance < np.nanpercentile(MinDistance[rgiExtrEval],75))/float(len(MinDistance)))*100.

    # calculate the AUC
    from sklearn.metrics import roc_auc_score
    testy=(testing_predictant >= np.sort(testing_predictant)[-extreme_nr])
    probs=(MinDistance-np.min(MinDistance)); probs=np.abs((probs/probs.max())-1)
    try:
        auc = roc_auc_score(testy, probs)
    except:
        auc=np.nan

    # Calculate the Average precision-recall score
    from sklearn.metrics import average_precision_score
    from sklearn import svm, datasets
    try:
        average_precision = average_precision_score(testy, probs)
    except:
        average_precision = np.nan

    # calculate precipitation pattern speciffic metrixc
    if len(testing_predictant) == len(rgrClustersFin[1]):
        # These skill scores are only calculate for "all day" weather typing
        if 'testing_predictant_GC' in locals():
    #         PR_Centroids_test = np.array([np.nanmean(testing_predictant_GC[rgrClustersFin[1] == ce, :], axis=0) for ce in np.unique(rgrClustersFin[1])])
    #         PR_anomaly = np.mean(np.abs(1-np.nanmean(PR_Centroids_test, axis=1)/np.nanmean(testing_predictant_GC)))
            PR_Centroids_test1 = np.array([np.nanmean(testing_predictant_GC[rgrClustersFin[1] == ce, :], axis=0) for ce in np.unique(rgrClustersFin[1])])
            PR_anomaly = np.nanmean(np.abs(1-PR_Centroids_test1/np.nanmean(testing_predictant_GC, axis=0)[None,:]))
    #         IntraClusterSTD = np.nanstd(np.nanmean(PR_Centroids_test,axis=1), axis=0)
    #         InterClusterSTD = np.array([np.std(np.nanmean(testing_predictant_GC[rgrClustersFin[1] == ce, :], axis=1), axis=0) for ce in np.unique(rgrClustersFin[1])])
    #         Intra_vs_Inter_Cluster_STD = IntraClusterSTD/np.mean(InterClusterSTD)
            IntraClusterSTD1 = np.nanstd(PR_Centroids_test1)
            InterClusterSTD1 = np.array([np.nanstd(testing_predictant_GC[rgrClustersFin[1] == ce, :]) for ce in np.unique(rgrClustersFin[1])])
            Intra_vs_Inter_Cluster_STD = IntraClusterSTD1/np.mean(InterClusterSTD1)
        else:
            PR_Centroids_test = np.array([np.nanmean(testing_predictant[rgrClustersFin[1] == ce]) for ce in np.unique(rgrClustersFin[1])])
            PR_anomaly = np.mean(np.abs(1-PR_Centroids_test/np.nanmean(testing_predictant)))
            IntraClusterSTD = np.nanstd(PR_Centroids_test, axis=0)
            InterClusterSTD = np.array([np.std(testing_predictant[rgrClustersFin[1] == ce], axis=0) for ce in np.unique(rgrClustersFin[1])])
            Intra_vs_Inter_Cluster_STD = IntraClusterSTD/np.mean(InterClusterSTD)
    else:
        PR_anomaly = None
        Intra_vs_Inter_Cluster_STD = None

    # print("--- Summary of performance ---")
    # print("    PSS: "+str(np.round(grPSS,2)))
    # print("    MRD: "+str(np.round(grMRD,2)))
    # print("    MRR: "+str(np.round(grMRR,2)))
    # print("    Excluded: "+str(np.round(grExluded,2)))
    # print("    AUC: "+str(np.round(auc,2)))
    # print("    APR: "+str(np.round(average_precision,2)))
    # print("------------------------------")

    XWT_output={'grClustersFin':rgrClustersFin, 
                'grEucledianDist':MinDistance, 
                'EucledianDistAllWTs':EucledianDist, 
                'grCorrelatio':MaxCorr,
                'grCorrelatioAllWTs':Correlation,
                'grPSS':grPSS,
                'grMRD':grMRD,
                'grMRR':grMRR,
                'APR':average_precision,
                'AUC':auc,
                'PEX':grExluded,
                'grExluded':grExluded,
                'PRanom':PR_anomaly,
                'InterVSIntra':Intra_vs_Inter_Cluster_STD}
    return XWT_output





def VarDecomp_4D(DATA,DimensionNames):
    # code is base on: https://link.springer.com/content/pdf/10.1007/s10584-006-9228-x.pdf
    
    #Main Terms
    DN=DimensionNames
    V=1./(DATA.shape[0])*np.nansum((np.nanmean(DATA, axis=(1,2,3))-np.nanmean(DATA))**2)
    E=1./(DATA.shape[1])*np.nansum((np.nanmean(DATA, axis=(0,2,3))-np.nanmean(DATA))**2)
    D=1./(DATA.shape[2])*np.nansum((np.nanmean(DATA, axis=(0,1,3))-np.nanmean(DATA))**2)
    S=1./(DATA.shape[3])*np.nansum((np.nanmean(DATA, axis=(0,1,2))-np.nanmean(DATA))**2)
    # Interaction terms
    VE=1./(DATA.shape[0]*DATA.shape[1])*np.nansum((np.nanmean(DATA, axis=(2,3))-np.nanmean(DATA, axis=(1,2,3))[:,None]-np.nanmean(DATA, axis=(0,2,3))[None,:]+np.nanmean(DATA))**2)
    VD=1./(DATA.shape[0]*DATA.shape[2])*np.nansum((np.nanmean(DATA, axis=(1,3))-np.nanmean(DATA, axis=(1,2,3))[:,None]-np.nanmean(DATA, axis=(0,1,3))[None,:]+np.nanmean(DATA))**2)
    VS=1./(DATA.shape[0]*DATA.shape[3])*np.nansum((np.nanmean(DATA, axis=(1,2))-np.nanmean(DATA, axis=(1,2,3))[:,None]-np.nanmean(DATA, axis=(0,1,2))[None,:]+np.nanmean(DATA))**2)
    ED=1./(DATA.shape[1]*DATA.shape[2])*np.nansum((np.nanmean(DATA, axis=(0,3))-np.nanmean(DATA, axis=(0,2,3))[:,None]-np.nanmean(DATA, axis=(0,1,3))[None,:]+np.nanmean(DATA))**2)
    ES=1./(DATA.shape[1]*DATA.shape[3])*np.nansum((np.nanmean(DATA, axis=(0,2))-np.nanmean(DATA, axis=(0,2,3))[:,None]-np.nanmean(DATA, axis=(0,1,2))[None,:]+np.nanmean(DATA))**2)
    DS=1./(DATA.shape[2]*DATA.shape[3])*np.nansum((np.nanmean(DATA, axis=(0,1))-np.nanmean(DATA, axis=(0,1,3))[:,None]-np.nanmean(DATA, axis=(0,1,2))[None,:]+np.nanmean(DATA))**2)

    VarianceComponents=(V,E,D,S,VE,VD,VS,ED,ES,DS)
    ComponentNames=(DN[0],DN[1],DN[2],DN[3],DN[0]+DN[1],DN[0]+DN[2],DN[0]+DN[3],DN[1]+DN[2],DN[1]+DN[3],DN[2]+DN[3])
    
    return VarianceComponents, ComponentNames




def Centroids_to_NetCDF(NetCDFname,
                        XWT_output,
                        LonWT,
                        LatWT,
                        DailyVarsTrain,
                        rgdTime,
                        TimeEval,
                        VarsJoint,
                        rgiExtremeDays):

    if os.path.exists(NetCDFname):
        os.remove(NetCDFname)
        
    XWT_output['grClustersFin'][0]
    LonWT
    LatWT
    DailyVarsTrain
    from datetime import datetime
    today = datetime.today()
    timeAct=np.array(range(len(rgdTime)))[(rgdTime.year >= TimeEval.year[0]) & (rgdTime.year <= TimeEval.year[-1])]
    
    CentroidsAct=XWT_output['grClustersFin'][0]
    CentroidsAct=np.reshape(CentroidsAct,(CentroidsAct.shape[0],LonWT.shape[0],LonWT.shape[1],DailyVarsTrain.shape[3]))
    CentroidsAct=np.moveaxis(CentroidsAct, 3, -3)
    
    # WRITE THE NETCDF FILE
    root_grp = Dataset(NetCDFname, 'w', format='NETCDF4')
    # dimensions
    root_grp.createDimension('time', None)
    root_grp.createDimension('rlon', LonWT.shape[1])
    root_grp.createDimension('rlat', LonWT.shape[0])
    root_grp.createDimension('XWTs', XWT_output['grClustersFin'][1].max()+1)
    root_grp.createDimension('vars', DailyVarsTrain.shape[3])
    root_grp.createDimension('svars', 4)
    
    # variables
    centroids = root_grp.createVariable('centroids', 'f4', ('XWTs','vars','rlat','rlon'),fill_value=-999999)
    euclideandistance = root_grp.createVariable('euclediandistance', 'f4', ('time','XWTs'),fill_value=-999999)
    rlat = root_grp.createVariable('rlat', 'f4', ('rlat','rlon',))
    rlon = root_grp.createVariable('rlon', 'f4', ('rlat','rlon',))
    time = root_grp.createVariable('time', 'f8', ('time',))
    
    
    # Variable Attributes
    centroids.units = ""
    centroids.long_name = 'XWT centroid'
    centroids.coordinates = "rlon rlat"
    centroids.cell_methods = "time: mean"
    
    euclideandistance.units = "[-]"
    euclideandistance.long_name = "Euclidean Distances"
    centroids.cell_methods = "time: mean"
    
    time.calendar = "gregorian"
    time.units = "days since "+str(TimeEval.year[0])+"-"+str(TimeEval.month[0])+"-1 12:00:00"
    time.standard_name = "time"
    time.long_name = "time"
    time.axis = "T"
    
    rlon.standard_name = "grid_longitude"
    rlon.units = "degrees"
    rlon.axis = "X"
    
    rlat.standard_name = "grid_latitude"
    rlat.units = "degrees"
    rlat.axis = "Y"
    
    root_grp.description = "XWT OUTPUT"
    root_grp.history = "Created " + today.strftime("%d/%m/%y")
    root_grp.XWT_variables = str(VarsJoint)
    root_grp.XWT_days=str(rgiExtremeDays)
    root_grp.XWT_days_index=str(XWT_output['grClustersFin'][1])
    root_grp.APR_score=str(XWT_output['APR'])
    root_grp.PEX_score=str(np.round(XWT_output['grExluded'],2))+' %'
    root_grp.contact = "prein@ucar.edu"
    
    # write data to netcdf
    rlat[:]=LatWT
    rlon[:]=LonWT
    centroids[:]=CentroidsAct
    euclideandistance[:]=XWT_output['EucledianDistAllWTs']
    time[:]=timeAct
    root_grp.close()
