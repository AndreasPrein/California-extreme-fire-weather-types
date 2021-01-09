#!/usr/bin/env python

# # NOAAXWTs_in_Model.ipynb

# 
# '''File name: XWTs_in_Model.ipynb
#     Author: Andreas Prein
#     E-mail: prein@ucar.edu
#     Date created: 04.05.2020
#     Date last modified: 04.05.2020
# 
#     ############################################################## 
#     Purpos:
# 
#     Load XWTs centroids that were created by running:
#     /glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/programs/XWTing/XWTs_search_ERA5.py
# 
#     1) Read the nescessary variables from the model in which XWTs should be searched for
#     2) Bring the centroids and the model data to the common coarser resolution grid
#     3) Preprocess the model data
#     3) Calculate Eucledian Distances for each day
#     4) Save the data
# 
# '''

# In[1]:


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
# import pickle
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
from scipy.signal import wiener, filtfilt, butter, gaussian, freqz
from scipy.ndimage import filters
import pickle
import time
import xarray as xr

from Functions_Extreme_WTs import XWT
from Functions_Extreme_WTs import MRR, MRD, perkins_skill
from Functions_Extreme_WTs import PreprocessWTdata
from Functions_Extreme_WTs import EucledianDistance

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


import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import numpy as np

def interp_weights(xy, uv,d=2):
    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts):
    return np.einsum('nj,nj->n', np.take(values, vtx), wts)


# ### Program Setup

RE = int(sys.argv[1])


# ---------
# Setup clustering algorithm
ClusterMeth='hdbscan'  # current options are ['HandK','hdbscan']
ClusterBreakup = 0     # breakes up clusters that are unproportionally large (only for hdbscan)
RelAnnom=1             # 1 - calculates daily relative anomalies
NormalizeData='D'      # normalize variables | options are  - 'C' - climatology
                                                        # - 'D' - daily (default)
                                                        # - 'N' - none
MinDistDD=7            # minimum nr of days between XWT events
RemoveAnnualCycl=1     # remove annual cycle in varaiables with 21 day moving average filter
# ---------
# In[2]:


Season='Annual' # ['AMJJAS', 'ONDJFM']

DW_Regions=['Bay_Area','Central_Coast','LA','Modoc','Northeast','San_Diego','Siearas_East','Siearas_West']
sSubregionPR='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/Shapefiles/' #+Regions[sRegion]

# SelVar = np.array([['RH2AVG','MF2AVG','V850'],
#          ['T2AVG','T2MIN','RH2AVG'],
#          ['MR2AVG','SLPAVG','WSPD10'],
#          ['T2MIN','U10AVG','U500'],
#          ['RH500','T2AVG','T500'],
#          ['MF2AVG','PWAVG','V200'],
#          ['T2MIN','T500','WSPD500'],
#          ['RH2AVG','V500','WSPD200']])

SelVar = np.array([['RH2AVG','MF2AVG','V850'],
         ['RH2AVG','SLPAVG','WSPD10'],
         ['MR2AVG','SLPAVG','WSPD10'],
         ['T2MIN','WSPD10','Z500'],
         ['RH2AVG','RH500','U200'],
         ['MF2AVG','PWAVG','V200'],
         ['RH500','T2MIN','V200'],
         ['RH2AVG','V500','WSPD200']])

E5Vars=['CAPE','CIN','LCL','LFC','MF2AVG','MF500','MF850','MR2AVG','MR500','MR850','PWAVG','RH2AVG',        'RH500','RH850','SLPAVG','T2AVG','T2MAX','T2MIN','T500','T850','U10AVG','U200','U500','U850',        'V10AVG','V200','V500','V850','WSPD10','WSPD200','WSPD500','WSPD850','Z500']
E20Vars = ['CAPE','CIN','LCL','LFC','MF2','MF500','FLX850','Q2','Q500','Q850','PW','RH2',           'RH500','RH850','PSL','T2','T2max','T2min','T500','T850','U10','U200','U500','U850',           'V10','V200','V500','V850','UV10','UV200','UV500','UV850','ZG500']


Nr_XWTs = [4,6,10,8,10,4,10,8]
smoothing_radius = 0.5

YYYY_stamp = '2001-2019'
Season = 'Annual'
FireObs = 'MODIS'
REANAL = 'ERA5'


# In[4]:
Season='Annual' # ['AMJJAS', 'ONDJFM']

Data_All={}
for re in [RE]: #range(len(DW_Regions)):
    Region=DW_Regions[re]
    print('    process '+Region)

    Variables=SelVar[re,:]
    ExtrNR = Nr_XWTs[re]
    iNrOfExtremes = ExtrNR
    s = "-"; VarString=s.join(Variables)
    VarsJoint = VarString
    
    XWTdata='/glade/campaign/mmm/c3we/prein/Papers/2019_Janice-CA-Fire-WTs/data/AUC-APR/'+'Clusters'+str(iNrOfExtremes)+'_'+Region+'_'+YYYY_stamp+'_'+VarsJoint+'_'+Season+'_'+\
    FireObs+'_'+REANAL+'_'+ClusterMeth+'_Breakup-'+str(ClusterBreakup)+'_RelAnnom-'+str(RelAnnom)+'_Norm-'+NormalizeData+\
    '_DDdist-'+str(MinDistDD)+'_RemAnnCyc-'+str(RemoveAnnualCycl)
    sDataDir='/glade/campaign/mmm/c3we/prein/Papers/2019_Janice-CA-Fire-WTs/data/XWT_20CR_v3/AUC-APR/'
    sClusterSave=sDataDir+'Clusters'+str(ExtrNR)+'_'+Region+'_'+YYYY_stamp+'_'+VarString+'_'+Season+'_'+FireObs+'_'+REANAL
    
    # ### Load the centroids and shape file
    print('    Restore: '+XWTdata)
    with open(XWTdata, 'rb') as handle:
        npzfile = pickle.load(handle)
    XWTclusters=npzfile['grClustersFin']['Full']
    XWTlat=npzfile['LatWT']
    XWTlon=npzfile['LonWT']
    XWTtime=npzfile['rgdTime']
    SpatialSmoothing=npzfile['SpatialSmoothing']
        
    print('        Load shapefile')
    sf = shp.Reader(sSubregionPR+Region)
    df = read_shapefile(sf)
    for sf in range(df.shape[0]):
        ctr = df['coords'][sf]
        if len(ctr) > 10000:
            ctr=np.array(ctr)[::100,:] # carsen the shapefile accuracy
        else:
            ctr=np.array(ctr)
        grPRregion=mplPath.Path(ctr)
    
    # ### Read ERA20C data
    ncid=Dataset('/gpfs/fs1/collections/rda/data/ds131.3/anl/anl_mean_2015_VGRD_pres.nc', mode='r')
    E20Lat=np.squeeze(ncid.variables['latitude'][:])
    E20Lon=np.squeeze(ncid.variables['longitude'][:])
    ncid.close()
    E20Lon[E20Lon>180]=E20Lon[E20Lon>180]-360
    E20Lon2D, E20Lat2D = np.meshgrid(E20Lon, E20Lat)
    
    rgrGridCells=[(E20Lon2D.ravel()[ii],E20Lat2D.ravel()[ii]) for ii in range(len(E20Lat2D.ravel()))]
    rgrSRactP=np.zeros((E20Lat2D.shape[0]*E20Lat2D.shape[1]))
    TMP=np.array(grPRregion.contains_points(rgrGridCells))
    rgrSRactP[TMP == 1]=1
    rgrSRactP=np.reshape(rgrSRactP, (E20Lat2D.shape[0], E20Lat2D.shape[1]))
    
    # get grid cells that cover the centroid
    N_XWT=np.max(XWTlat)
    S_XWT=np.min(XWTlat)
    E_XWT=np.max(XWTlon)
    W_XWT=np.min(XWTlon)
    
    # ERA-20C area that contains this region
    E20_W=np.where(np.max((E20Lon2D <= W_XWT), axis=0) == 1)[0][-1]+1 #-1
    E20_E=np.where(np.max((E20Lon2D <= E_XWT), axis=0) == 1)[0][-1]+1
    E20_S=np.where(np.max((E20Lat2D <= S_XWT), axis=1) == 1)[0][0]-1
    E20_N=np.where(np.max((E20Lat2D <= N_XWT), axis=1) == 1)[0][0] #-1
    
    # # Cut out Lat and Lon
    E20Lon2D=E20Lon2D[E20_N:E20_S,E20_W:E20_E]
    E20Lat2D=E20Lat2D[E20_N:E20_S,E20_W:E20_E]
    
    print('    process the regridding weights')
    [Y,X]=(XWTlon, XWTlat)
    [Yi,Xi]=(E20Lon2D,E20Lat2D)
    xy=np.zeros([X.shape[0]*X.shape[1],2])
    xy[:,0]=Y.flatten()
    xy[:,1]=X.flatten()
    uv=np.zeros([Xi.shape[0]*Xi.shape[1],2])
    uv[:,0]=Yi.flatten()
    uv[:,1]=Xi.flatten()
    #Computed once and for all !
    vtx, wts = interp_weights(xy, uv)
    
    # Read E20C
    dStartDayE20=datetime.datetime(1900, 1, 1,0)
    dStopDayE20=datetime.datetime(2015, 12, 31,23)
    E20_time = pd.date_range(dStartDayE20, end=dStopDayE20, freq='d')
    E20YYYY=np.unique(E20_time.year)
    
    # ### Bring data to same grid
    if XWTlat.shape[0]*XWTlon.shape[0] > E20Lat2D.shape[0]*E20Lat2D.shape[1]:
        # we remap the centroids to the model grid
        Clusters=np.reshape(XWTclusters[0],(XWTclusters[0].shape[0],XWTlat.shape[0],XWTlon.shape[1],len(Variables)))
        RemappedClusters=np.zeros((Clusters.shape[0],E20Lat2D.shape[0],E20Lat2D.shape[1],Clusters.shape[3]))

        for cl in range(Clusters.shape[0]):
            for va in range(len(Variables)):
                values=Clusters[cl,:,:,va]
                valuesi=interpolate(values.flatten(), vtx, wts)
                RemappedClusters[cl,:,:,va]=valuesi.reshape(Xi.shape[0],Xi.shape[1])
        # we have to renormalize the clusters
        Normalized=(RemappedClusters-np.mean(RemappedClusters, axis=(1,2))[:,None,None,:])
        Normalized=Normalized/np.std(Normalized, axis=(1,2))[:,None,None,:]
        ClusterPrepared=(np.reshape(Normalized, (Normalized.shape[0],Normalized.shape[1]*Normalized.shape[2]*Normalized.shape[3])),XWTclusters[1])
    else:
        print('    Remapping the model data to the centroid grid is not implemented yet!')
        stop()
    for me in range(80):
        start = time.process_time()
        print('        member '+str(me+1).zfill(3))
        E20_DATA=np.zeros((len(E20_time),E20Lon2D.shape[0],E20Lon2D.shape[1], len(Variables))); E20_DATA[:]=np.nan
        for yy in range(len(E20YYYY)):
            DDact=(E20_time.year == E20YYYY[yy])
            for va in range(len(Variables)):
                VARact = E20Vars[E5Vars.index(Variables[va])]
                FILE='/glade/campaign/mmm/c3we/prein/NOAA-20C/daymean/'+VARact+'/'+VARact+'.'+str(E20YYYY[yy])+'_mem'+str(me+1).zfill(3)+'.nc'
#                 ds = xr.open_dataset(FILE)
#                 E20_DATA[DDact,:,:,va]=ds[VARact].sel(rlon=slice(E20_W, E20_E), rlat=slice(E20_N, E20_S)).values
                ncid=Dataset(FILE, mode='r') # open the netcdf file
                try:
                    E20_DATA[DDact,:,:,va]=np.squeeze(ncid.variables[VARact][:,E20_N:E20_S,E20_W:E20_E])
                except:
                    DATA_act = np.squeeze(ncid.variables[VARact][:,E20_N:E20_S,E20_W:E20_E])
                    E20_DATA[np.where(DDact == True)[0][:DATA_act.shape[0]],:,:,va]=DATA_act
                ncid.close()
        
        # your code here    
        print(time.process_time() - start)

        DailyVarsEvalNorm=PreprocessWTdata(E20_DATA,                      # WT data [time,lat,lon,var]
                               RelAnnom=RelAnnom,                         # calculate relative anomalies [1-yes; 0-no]
                               SmoothSigma=smoothing_radius,              # Smoothing stddev (Gaussian smoothing)
                               RemoveAnnualCycl=RemoveAnnualCycl,         # remove annual cycle [1-yes; 0-no]
                               NormalizeData=NormalizeData)               # normalize data [1-yes; 0-no]

        EucledianDist, Correlation =EucledianDistance(DailyVarsEvalNorm,
                                                      ClusterPrepared)
        EucledianDist_orig=np.copy(EucledianDist)
        EucledianDist=np.min(EucledianDist,axis=1)
        
        EDspace=np.linspace(EucledianDist.min(),np.percentile(EucledianDist,0.2),100)
        ChangeInFr=np.zeros((len(EDspace)))
        RP=np.zeros((len(EDspace)))
        ClimLen=40

        for ed in range(len(EDspace)):
            TEST=(EucledianDist < EDspace[ed])
            ExtrFrY=float(np.sum(TEST))/len(E20YYYY)
            ExpectedFrequency=np.cumsum(np.array([ExtrFrY]*len(E20YYYY)))
            RP[ed]=ExtrFrY
            TEST=np.array([np.sum(TEST[E20_time.year == E20YYYY[yy]]) for yy in range(len(E20YYYY))])
            TEST=np.cumsum(TEST, dtype=float)

            Ref=(TEST[ClimLen]-TEST[0])
            Cur=(TEST[-1]-TEST[-ClimLen])
            ChangeInFr[ed]=((Cur-Ref)/float(Ref))*100

        # In[37]:
        print('    calculate random chance')
        try: RandomChance
        except:
            BSnr=1000
            RandomChance=np.zeros((len(RP),BSnr))
            for pe in range(len(RP)):
                NrOfEvents=(RP*len(E20YYYY))[pe]
                for bs in range(BSnr):
                    RandomDays=np.sort(np.array([random.randint(0, len(E20YYYY)-1) for ii in range(int(NrOfEvents))]))
                    RandRecord=np.zeros((len(E20YYYY)))
                    for ii in range(len(RandomDays)):
                        RandRecord[RandomDays[ii]]=RandRecord[RandomDays[ii]]+1
                    Ref=np.sum(RandRecord[:ClimLen])
                    Cur=np.sum(RandRecord[-ClimLen:])
                    RandomChance[pe,bs]=((Cur-Ref)/float(Ref))*100

        Data_reg={'EucledianDist_orig':EucledianDist_orig,
                  'ChangeInFr':ChangeInFr,
                  'RP':RP,
                  'RandomChance':RandomChance,
                  'E20_time':E20_time}
        Data_All[Region+'_'+str(me).zfill(3)]=Data_reg

# Save the data
afile = open(sDataDir+'20CR_v3_XWTs_Eucledian-Distances_'+DW_Regions[RE]+'.pkl', 'wb')
pickle.dump(Data_All, afile)
afile.close()



