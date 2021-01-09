#!/usr/bin/env python
'''File name: ExtremeEvent-WeatherTyping.py
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 16.04.2018
    Date last modified: 16.04.2018

    ############################################################## 
    Purpos:
    Classifies weather patterns that cause precipitation extremes

    1) read in shape file for area under cosideration

    2) read in precipitation data from PRISM

    3) identify the N-days that had highest rainfall records

    4) read in ERA-Interim data for these days

    5) remove the 30-year mooving average from time series and
       normalize the variables

    5) run clustering algorithm on extreme WT patterns

    6) search for the extreme WT centroids in the full record


'''

import matplotlib
gui_env = ['TKAgg','GTKAgg','Qt4Agg','WXAgg']
for gui in gui_env:
    try:
        print( "testing", gui)
        matplotlib.use(gui,warn=False, force=True)
        from matplotlib import pyplot as plt
        break
    except:
        continue
print( "Using:",matplotlib.get_backend())


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
import shapefile as shp
import sys
from itertools import combinations 

from Functions_Extreme_WTs import XWT
from Functions_Extreme_WTs import MRR, MRD, perkins_skill

import warnings
warnings.filterwarnings("ignore")

# ###################################################
# This information comes from the setup file

from XWTs_search_ERA5 import rgdTime, iMonths, sPlotDir, sDataDir, sSubregionPR, rgsWTvars, VarsFullName,rgsWTfolders, rgrNrOfExtremes, WT_Domains, DomDegreeAdd, Annual_Cycle, SpatialSmoothing, Metrics, Dimensions, DW_Regions, FireObs, ClusterMeth, ClusterBreakup, RelAnnom, NormalizeData, MinDistDD, RemoveAnnualCycl
# from XWTs_search_ERA5_preselect import rgdTime, iMonths, sPlotDir, sDataDir, sSubregionPR, rgsWTvars, VarsFullName,rgsWTfolders, rgrNrOfExtremes, WT_Domains, DomDegreeAdd, Annual_Cycle, SpatialSmoothing, Metrics, Dimensions, DW_Regions, FireObs

rgsWTfolders=[rgsWTfolders[0]+str(rgsWTvars[va])+'/' for va in range(len(rgsWTvars))]

# create all possible combinations of variables
iVariables=range(len(VarsFullName))
Combinations1=np.array(list(combinations(iVariables, 1)))
Combinations2=np.squeeze(np.array(list(combinations(iVariables, 2))))
Combinations3=np.squeeze(np.array(list(combinations(iVariables, 3))))
Combinations4=np.squeeze(np.array(list(combinations(iVariables, 4))))
Combinations=list(Combinations1)+list(Combinations2)+list(Combinations3) #+list(Combinations4)
print('Number of var. combinations '+str(len(Combinations)))

print('Number of arguments:', len(sys.argv), 'arguments.')
iRegion=int(sys.argv[1])
sRegion=DW_Regions[iRegion]
sSubregionPR=sSubregionPR+sRegion
print('---- Process Region '+sRegion+' ----')

# create nessesary directories
if not os.path.exists(sDataDir):
    os.makedirs(sDataDir)
if not os.path.exists(sPlotDir):
    os.makedirs(sPlotDir)
sRegion=sRegion.replace('/','-')

ss='-'
sMonths=ss.join([str(iMonths[ii]) for ii in range(len(iMonths))])
print('    Months: '+sMonths)

# ###################################################
# use setup to generate data
rgiYears=np.unique(rgdTime.year)
rgdTimeMM = pd.date_range(datetime.datetime(rgdTime.year[0], rgdTime.month[0], rgdTime.day[0],0), end=datetime.datetime(rgdTime.year[-1], rgdTime.month[-1], rgdTime.day[-1],0), freq='m')
AllMonths=rgdTime
YYYY_stamp=str(rgdTime.year[0])+'-'+str(rgdTime.year[-1])
rgiSeasonWT=np.isin(rgdTime.month, iMonths)
rgdTime=rgdTime[rgiSeasonWT]

SPLIT=np.where(rgdTime.year < rgiYears[int(np.round(len(rgiYears)/2))])[0][-1]
    
SkillScores_All=np.zeros((len(Combinations),len(rgrNrOfExtremes),len(WT_Domains),len(Annual_Cycle),len(SpatialSmoothing),2,len(Metrics))); SkillScores_All[:]=np.nan

### CHECK IF DATA IS ALREADY PROCESSED ###
SaveStats=sDataDir+sRegion+'_'+YYYY_stamp+'-'+sMonths+'_'+FireObs+'_ERA5.npz'
if os.path.isfile(SaveStats) == 0:
    # ###################################################
    print('    Read the daily Fire data')
    
    
    ncid=Dataset('/glade/campaign/mmm/c3we/prein/observations/Fire-observations/Fire_GWIS/data_20020702/GriddedData/MODIS_BA_GLOBAL_1_2016_01_gridded.nc', mode='r') # open the netcdf file
    rgrLatPR=np.squeeze(ncid.variables['rlat'][:])
    rgrLonPR=np.squeeze(ncid.variables['rlon'][:])
    ncid.close()
    rgrGridCells=[(rgrLonPR.ravel()[ii],rgrLatPR.ravel()[ii]) for ii in range(len(rgrLonPR.ravel()))]
    rgrSRactP=np.zeros((rgrLonPR.shape[0]*rgrLonPR.shape[1]))
    sf = shp.Reader(sSubregionPR)
    from Functions_Extreme_WTs import read_shapefile
    df = read_shapefile(sf)
    for sf in range(df.shape[0]):
        ctr = df['coords'][sf]
        if len(ctr) > 10000:
            ctr=np.array(ctr)[::100,:] # carsen the shapefile accuracy
        else:
            ctr=np.array(ctr)
        grPRregion=mplPath.Path(ctr)
        TMP=np.array(grPRregion.contains_points(rgrGridCells))
        rgrSRactP[TMP == 1]=1

    rgrSRactP=np.reshape(rgrSRactP, (rgrLatPR.shape[0], rgrLatPR.shape[1]))
    rgiSrPR=np.array(np.where(rgrSRactP == True))
    iLatMaxP=rgiSrPR[0,:].max()+1
    iLatMinP=rgiSrPR[0,:].min()
    iLonMaxP=rgiSrPR[1,:].max()+1
    iLonMinP=rgiSrPR[1,:].min()
    rgrPRdata=np.zeros((sum(rgiSeasonWT),iLatMaxP-iLatMinP,iLonMaxP-iLonMinP))
        
    if FireObs == 'MODIS':
        for mo in range(len(rgdTimeMM)):
            rgiDD=np.where(((rgdTimeMM.year[mo] == rgdTime.year) & (rgdTime.month == rgdTimeMM.month[mo]) & np.isin(rgdTimeMM.month[mo], iMonths)))[0]
            if len(rgiDD) != 0:
                ncid=Dataset('/glade/campaign/mmm/c3we/prein/observations/Fire-observations/Fire_GWIS/data_20020702/GriddedData/MODIS_BA_GLOBAL_1_'+str(rgdTimeMM.year[mo])+'_'+str("%02d" % rgdTimeMM.month[mo])+'_gridded.nc', mode='r')
                rgrPRdata[rgiDD,:,:]=np.squeeze(ncid.variables['BurnedArea'][:,iLatMinP:iLatMaxP,iLonMinP:iLonMaxP])
                ncid.close()
        
        rgrPRdata[rgrPRdata<0] = np.nan
        rgiSRgridcells=rgrSRactP[iLatMinP:iLatMaxP,iLonMinP:iLonMaxP].astype('int')
        rgrPRdata=np.nanmean(rgrPRdata[:,(rgiSRgridcells==1)], axis=(1))

    elif FireObs == 'Parks':
        SheanDataDir='/glade/campaign/mmm/c3we/prein/observations/Fire-observations/Sean_Parks/SubRegion_DailyBurnedArea/'
        FILE=SheanDataDir+DW_Regions[iRegion]+'.txt'
        rgrPRrecords=np.array(pd.read_csv(FILE, delim_whitespace=True, header = None))[:,0]
        
        # rgiExtrTrain=ExtremeDays(rgrPRrecords,iNrOfExtremes,MinDistDD)
        # rgiExtremeDays=rgdTime[rgiExtrTrain]
        rgrPRdata=np.copy(rgrPRrecords)
        rgrPRdata[rgrPRdata<0] = np.nan

    # We read in ERA5 data for the largest region and cut it to fit smaller regions
    DomDelta=np.max(DomDegreeAdd)
    Wlon=ctr[:,0].min()
    Elon=ctr[:,0].max()
    Nlat=ctr[:,1].max()
    Slat=ctr[:,1].min()
    DomainWT=np.array([[Elon+DomDelta,Slat-DomDelta],
                       [Wlon-DomDelta,Slat-DomDelta],
                       [Wlon-DomDelta,Nlat+DomDelta],
                       [Elon+DomDelta,Nlat+DomDelta],
                       [Elon+DomDelta,Slat-DomDelta]])
    grWTregion=mplPath.Path(DomainWT)

    # ###################################################
    #         Read the ERA-Interim grid and data
    print('    Read the ERA5 data specific for the region')
    from Functions_Extreme_WTs import ReadERA5
    DailyVarsLargeDom=ReadERA5(grWTregion,      # shapefile with WTing region
                       rgdTime,         # time period for WTing
                       iMonths,         # list of months that should be considered
                       rgsWTfolders,    # directories containing WT files
                       rgsWTvars)       # netcdf variable names of WT variables

    # for CAPE and CIN replace NaN with 0
    VARSact=['CAPE','CIN']
    for vv in range(len(VARSact)):
        try:
            IND=VarsFullName.index(VARSact[vv])
            TMP=DailyVarsLargeDom[0][:,:,:,IND]; TMP[np.isnan(TMP)]=0
            DailyVarsLargeDom[0][:,:,:,IND]=TMP
        except:
            pass
    # ###################################################
    RestartFolder = sDataDir+'/restarts/'
    if os.path.isdir(RestartFolder) != 1:
        subprocess.call(["mkdir","-p",RestartFolder])
    RestartFile = RestartFolder + 'Restartfile_'+sRegion+'_'+YYYY_stamp+'-'+sMonths+'_'+FireObs+'_ERA5_'+ClusterMeth+'_Breakup-'+str(ClusterBreakup)+'_RelAnnom-'+str(RelAnnom)+'_Norm-'+NormalizeData+'_DDdist-'+str(MinDistDD)+'_RemAnnCyc-'+str(RemoveAnnualCycl)+'.npz'
    if os.path.isfile(RestartFile) == 0:
        re0 = 0
        ss0 = 0
        ne0 = 0
        sm0 = 0
        ac0 = 0
    else:
        print('--------------')
        print('Restart calculation from file '+RestartFile)
        print('--------------')
        DATA = np.load(RestartFile)
        re0 = DATA['re0']
        ss0 = DATA['ss0']
        ne0 = DATA['ne0']
        sm0 = DATA['sm0']
        ac0 = DATA['ac0']
        SkillScores_All = DATA['SkillScores_All']

    for re in range(re0,len(WT_Domains),1):
        print( '    ------')
        print( '    Domain '+WT_Domains[re])
        DeltaX=np.max(DomDegreeAdd)-DomDegreeAdd[re]
        if DeltaX != 0:
            DomainWT=np.array([[Elon+DomDegreeAdd[re],Slat-DomDegreeAdd[re]],
                       [Wlon-DomDegreeAdd[re],Slat-DomDegreeAdd[re]],
                       [Wlon-DomDegreeAdd[re],Nlat+DomDegreeAdd[re]],
                       [Elon+DomDegreeAdd[re],Nlat+DomDegreeAdd[re]],
                       [Elon+DomDegreeAdd[re],Slat-DomDegreeAdd[re]]])

            grWTregion=mplPath.Path(DomainWT)
            rgrGridCells=[(DailyVarsLargeDom[1].ravel()[ii],DailyVarsLargeDom[2].ravel()[ii]) for ii in range(len(DailyVarsLargeDom[1].ravel()))]
            rgrSRact=np.array(grWTregion.contains_points(rgrGridCells)); rgrSRact=np.reshape(rgrSRact, (DailyVarsLargeDom[1].shape[0], DailyVarsLargeDom[1].shape[1]))
            rgiSrWT=np.array(np.where(rgrSRact == True))
            iLatMax=rgiSrWT[0,:].max()
            iLatMin=rgiSrWT[0,:].min()
            iLonMax=rgiSrWT[1,:].max()
            iLonMin=rgiSrWT[1,:].min()

            DailyVars=DailyVarsLargeDom[0][:,iLatMin:iLatMax,iLonMin:iLonMax,:]
        else:
            DailyVars=DailyVarsLargeDom[0]
        # perform split sample statistic
        for ss in range(ss0,2,1):
            print('    Split Sample Nr. '+str(ss+1))
            if ss == 0:
                DailyVarsTrain=DailyVars[:SPLIT,:]
                DailyVarsEval=DailyVars[-SPLIT:,:]
                Ptrain=rgrPRdata[:SPLIT]
                Peval=rgrPRdata[-SPLIT:]
                TimeTrain=rgdTime[:SPLIT]
                TimeEval=rgdTime[-SPLIT:]
            else:
                DailyVarsTrain=DailyVars[-SPLIT:,:]
                DailyVarsEval=DailyVars[:SPLIT,:]
                Ptrain=rgrPRdata[-SPLIT:]
                Peval=rgrPRdata[:SPLIT]
                TimeTrain=rgdTime[-SPLIT:]
                TimeEval=rgdTime[:SPLIT]

            for ne in range(ne0,len(rgrNrOfExtremes),1):
                DailyVarsAct=np.copy(DailyVarsTrain)
                print( '        '+str(rgrNrOfExtremes[ne])+' EXTREMES')
                iNrOfExtremes=rgrNrOfExtremes[ne]   # we consider the N highest rainfall extremes
                rgiSRgridcells=rgrSRactP[iLatMinP:iLatMaxP,iLonMinP:iLonMaxP].astype('int')
                rgrPRrecords=Ptrain
                rgrPReval=Peval
    
                # Test effect of spatial smoothing
                for sm in range(sm0,len(SpatialSmoothing),1):
                    # annual cycle treatment
                    for ac in range(ac0,len(Annual_Cycle),1):
                        print( '            Loop over variable permutations')
                        for va1 in range(len(Combinations)):
#                             print('    process var '+str(va1)+' of '+str(len(Combinations)))
                            XWT_output=XWT(DailyVarsTrain[:,:,:,Combinations[va1]],
                                           DailyVarsEval[:,:,:,Combinations[va1]],
                                           rgrPRrecords,
                                           rgrPReval,
                                           TimeTrain,
                                           TimeEval,
                                           rgrNrOfExtremes[ne],
                                           SpatialSmoothing[sm],
                                           ClusterMeth=ClusterMeth,
                                           ClusterBreakup=ClusterBreakup,
                                           RelAnnom=RelAnnom,
                                           NormalizeData=NormalizeData,
                                           MinDistDD=MinDistDD,
                                           RemoveAnnualCycl=RemoveAnnualCycl)                            
                            
                            
                            if XWT_output != None:
                                SkillScores_All[va1, ne, re, ac, sm, ss, Metrics.index('PSS')]=XWT_output['grPSS'] # Perkins Skill Score
                                SkillScores_All[va1, ne, re, ac, sm, ss, Metrics.index('MRD')]=XWT_output['grMRD'] # Mean relative difference
                                SkillScores_All[va1, ne, re, ac, sm, ss, Metrics.index('MRR')]=XWT_output['grMRR'] # Mean Rank Ratio
                                SkillScores_All[va1, ne, re, ac, sm, ss, Metrics.index('APR')]=XWT_output['APR'] # Average precision-recall score
                                SkillScores_All[va1, ne, re, ac, sm, ss, Metrics.index('PEX')]=XWT_output['PEX'] # Percent of points excluded for ED larger than the 75 percentile
                                SkillScores_All[va1, ne, re, ac, sm, ss, Metrics.index('AUC')]=XWT_output['AUC'] # Area under the ROC curve
                        print('    save restart file')
                        DATA = np.savez(RestartFile,
                                        re0 = re,
                                        ss0 = ss,
                                        ne0 = ne,
                                        sm0 = sm,
                                        ac0 = ac,
                                        SkillScores_All = SkillScores_All)
                print( ' ')

#     plt.scatter(SkillScores_All[:,:,:,:,:,:,3],SkillScores_All[:,:,:,:,:,:,5]); plt.show()
    np.savez(SaveStats,
             SkillScores_All=SkillScores_All, 
             Combinations=Combinations, 
             rgsWTvars=VarsFullName,
             rgrNrOfExtremes=rgrNrOfExtremes,
             WT_Domains=WT_Domains,
             Annual_Cycle=Annual_Cycle,
             SpatialSmoothing=SpatialSmoothing,
             Metrics=Metrics,
             Dimensions=Dimensions)
else:
    print('    Load: '+SaveStats)
    DATA=np.load(SaveStats)
    SkillScores_All=DATA['SkillScores_All']
    Combinations=DATA['Combinations']
    VarsFullName=DATA['rgsWTvars']
    rgrNrOfExtremes=DATA['rgrNrOfExtremes']
    WT_Domains=DATA['WT_Domains']
    Annual_Cycle=DATA['Annual_Cycle']
    SpatialSmoothing=DATA['SpatialSmoothing']
    Metrics=DATA['Metrics']
    Dimensions=DATA['Dimensions']

# Find optimum and print best setting
Metrics=list(Metrics)
Scores=[Metrics.index('APR'),Metrics.index('AUC')] #[Metrics.index('APR'),Metrics.index('AUC')]
Mean_SS=np.nanmean(SkillScores_All[:,:,:,:,:,:,Scores], axis=(5,6))
iOpt=np.where(np.nanmax(Mean_SS) == Mean_SS)


print(' ')
print('====================================')
print('======    OPTIMAL SETTINGS    ======')
print('VARIABLES')
for va in range(len(Combinations[iOpt[0][0]])):
    print('    '+VarsFullName[int(Combinations[iOpt[0][0]][va])])
print('Extreme Nr     : '+str(rgrNrOfExtremes[iOpt[1][0]]))
print('Domain Size    : '+str(WT_Domains[iOpt[2][0]]))
print('Annual Cy. Rem.: '+str(Annual_Cycle[iOpt[3][0]]))
print('Smoothing      : '+str(SpatialSmoothing[iOpt[4][0]]))
print('APR            : '+str(np.round(SkillScores_All[iOpt[0][0],iOpt[1][0],iOpt[2][0],iOpt[3][0],iOpt[4][0],:,Metrics.index('APR')],2)))
print('AUC            : '+str(np.round(SkillScores_All[iOpt[0][0],iOpt[1][0],iOpt[2][0],iOpt[3][0],iOpt[4][0],:,Metrics.index('AUC')],2)))
print('Average Score  : '+str(np.round(np.nanmax(Mean_SS),2))+' | 1 is perfect')
print('====================================')

stop()
# PlotFile=sRegion+'_XWT_Search-Optimum_'+YYYY_stamp+'_'+sMonths+'.pdf'
# from Functions_Extreme_WTs import SearchOptimum_XWT
# SearchOptimum_XWT(PlotFile,
#                  sPlotDir,
#                  SkillScores_All,
#                  GlobalMinimum1,
#                  GlobalMinimum2,
#                  Optimum,
#                  VariableIndices,
#                  Dimensions,
#                  Metrics,
#                  VarsFullName,
#                  ss,
#                  rgrNrOfExtremes,
#                  WT_Domains,
#                  Annual_Cycle,
#                  SpatialSmoothing)


