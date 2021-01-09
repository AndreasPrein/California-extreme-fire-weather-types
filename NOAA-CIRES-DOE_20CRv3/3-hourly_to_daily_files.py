#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python


# # 3-hourly_to_daily_files

# In[3]:


'''
    File name: 3-hourly_to_daily_files
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 23.07.2020
    Date last modified: 23.07.2020

    ##############################################################
    Purpos:

    1) Reads in 3-hourly original data from the NOAA-CIRES-DOE 20CRv3
       The original data is here - https://portal.nersc.gov/archive/home/projects/incite11/www/20C_Reanalysis_version_3/everymember_anal_netcdf/subdaily
       Chi-Fan Shih <chifan@ucar.edu> helped transfering the data via Globus

    2) untarr the data and calculate daily average NetCDF files for each of the 80 members in annual files


'''


# In[4]:


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
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from pylab import *
import string
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import shapefile
import shapely.geometry
# import descartes
import shapefile
import math
from scipy.stats.kde import gaussian_kde
from math import radians, cos, sin, asin, sqrt
from scipy import spatial
import scipy.ndimage
import matplotlib.path as mplPath
from scipy.interpolate import interp1d
import time
from math import atan2, degrees, pi
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
# import SkewT
import csv
import pygrib
from scipy import interpolate
import re

# In[21]:


# ________________________________________________________________________
# ________________________________________________________________________
#             USER INPUT SECTION

DataDir='/glade/campaign/mmm/c3we/prein/NOAA-20C/original_data/'
SaveDir='/glade/campaign/mmm/c3we/prein/NOAA-20C/daymean/'

VAR_SEL = str(sys.argv[1])

ConstantVars='/gpfs/fs1/collections/rda/data/ds131.3/anl/anl_mean_2015_VGRD_pres.nc'
Labels = ['OrigVar','TarVar','daystat','longname','unit']
VARS=[['PRES','PRES','mean','surface pressure','hPa'],
      ['PRMSL','PSL','mean','sea level pressure','hPa'],
      ['PWAT','PW','mean','precipitable water','mm'],
      ['SPFH2m','Q2','mean','2m specific humidity', 'kg kg-1'],
      ['SPFH500','Q500','mean','500 hPa specific humidity','kg kg-1'],
      ['TMP2m','T2','mean','2m mean air temperature', 'K'],
      ['TMP2m','T2min','min','2m minimum air temperature', 'K'],
      ['TMP2m','T2max','max','2m maximum air temperature', 'K'],
      ['TMP500','T500','mean','500 hPa mean air temperature', 'K'],
      ['UGRD10m','U10','mean','10m zonal wind speed', 'm s-1'],
      ['UGRD200','U200','mean','200hPa zonal wind speed', 'm s-1'],
      ['UGRD500','U500','mean','500hPa zonal wind speed', 'm s-1'],
      ['UGRD850','U850','mean','850hPa zonal wind speed', 'm s-1'],
      ['VGRD10m','V10','mean','10m meridional wind speed', 'm s-1'],
      ['VGRD200','V200','mean','200hPa meridional wind speed', 'm s-1'],
      ['VGRD500','V500','mean','500hPa meridional wind speed', 'm s-1'],
      ['VGRD850','V850','mean','850hPa meridional wind speed', 'm s-1'],
      ['HGT500','ZG500','mean','500hPa geopotential height', 'm']]
VARS_In = [VARS[ii][1] for ii in range(len(VARS))]

iVAR = VARS_In.index(VAR_SEL)

rgdTime = pd.date_range(datetime.datetime(1900, 1, 1,0), end=datetime.datetime(2015, 12, 31,23), freq='d')
YYYYall = np.unique(rgdTime.year)


# In[10]:


# ________________________________________________________________________
# ________________________________________________________________________
#              READ IN THE COORDINATES
ncid=Dataset(ConstantVars, mode='r')
rgrLat=np.squeeze(ncid.variables['latitude'][:])
rgrLon=np.squeeze(ncid.variables['longitude'][:])
ncid.close()


# In[ ]:

for va in [iVAR]: #range(len(VARS)):
    print('    processing '+VARS[va][0])
    OutDir = SaveDir+VARS[va][1]+'/'
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    
    for yy in range(len(YYYYall)):
        print('        '+str(YYYYall[yy]))
        # unpack the data
        FILES = glob.glob(OutDir+VARS[va][1]+'.'+str(YYYYall[yy])+"*.nc")
        if len(FILES) < 80:
            pp = subprocess.Popen("tar -xf "+DataDir+VARS[va][0]+'_'+str(YYYYall[yy])+'.tar -C '+OutDir, shell=True)
            pp.wait()
            OutDirEns = OutDir+str(YYYYall[yy])
            for em in range(80):
                FileACT = OutDirEns+'/'+VARS[va][0]+'.'+str(YYYYall[yy])+'_mem'+str(em+1).zfill(3)+'.nc'
                File_fin = OutDir+VARS[va][1]+'.'+str(YYYYall[yy])+'_mem'+str(em+1).zfill(3)+'.nc'
                File=File_fin+'_copy'
                if os.path.exists(File_fin) == 0:
                    match = re.match(r"([a-z]+)([0-9]+)", VARS[va][0], re.I)
                    if match:
                        items = match.groups()
                        NC_Name = items[0]
                        if (VARS[va][0] == 'TMP2m') | (VARS[va][0] == 'UGRD10m') | (VARS[va][0] == 'VGRD10m'):
                            NC_Name = VARS[va][0]
                    else:
                        NC_Name = VARS[va][0]
                    try:
                        ncid=Dataset(FileACT, mode='r')
                        rgrDATA=np.squeeze(ncid.variables[NC_Name][:])
                        time_act=np.squeeze(ncid.variables['time'][:])[4::8]
                        ncid.close()
                    except:
                        ncid=Dataset(OutDirEns+'/'+VARS[va][0]+'.'+str(YYYYall[yy])+'_mem'+str(em).zfill(3)+'.nc', mode='r')
                        rgrDATA=np.squeeze(ncid.variables[NC_Name][:])
                        time_act=np.squeeze(ncid.variables['time'][:])[4::8]
                        ncid.close()
                    

                    if VARS[va][2] == 'mean':
                        rgrDATA = np.mean(np.reshape(rgrDATA, (int(rgrDATA.shape[0]/8),8,rgrDATA.shape[1],rgrDATA.shape[2])), axis=1)
                    elif VARS[va][2] == 'min':
                        rgrDATA = np.min(np.reshape(rgrDATA, (int(rgrDATA.shape[0]/8),8,rgrDATA.shape[1],rgrDATA.shape[2])), axis=1)
                    elif VARS[va][2] == 'max':
                        rgrDATA = np.max(np.reshape(rgrDATA, (int(rgrDATA.shape[0]/8),8,rgrDATA.shape[1],rgrDATA.shape[2])), axis=1)
    #                 except:
    #                     rgrDATA = rgrDATA[-int(rgrDATA.shape[0]/8)*8:,:]
    #                     rgrDATA = np.mean(np.reshape(rgrDATA, (int(rgrDATA.shape[0]/8),8,rgrDATA.shape[1],rgrDATA.shape[2])), axis=1)

                    # ________________________________________________________________________
                    # write the netcdf
#                     print( '        ----------------------')
#                     print( '        Save data to '+File_fin)
                    root_grp = Dataset(File, 'w', format='NETCDF4')
                    # dimensions
                    root_grp.createDimension('time', None)
                    root_grp.createDimension('rlon', rgrLon.shape[0])
                    root_grp.createDimension('rlat', rgrLat.shape[0])
                    # variables
                    lat = root_grp.createVariable('lat', 'f4', ('rlat',))
                    lon = root_grp.createVariable('lon', 'f4', ('rlon',))
                    time = root_grp.createVariable('time', 'f8', ('time',))
                    UV = root_grp.createVariable(VARS[va][1], 'f4', ('time','rlat','rlon',),fill_value=-99999)

                    time.calendar = "proleptic_gregorian"
                    time.units = "hours since 1800-01-01 00:00:00.0"
                    time.standard_name = "time"
                    time.long_name = "time"

                    lon.standard_name = "longitude"
                    lon.long_name = "longitude"
                    lon.units = "degrees_east"

                    lat.standard_name = "latitude"
                    lat.long_name = "latitude"
                    lat.units = "degrees_north"

                    UV.standard_name = VARS[va][3]
                    UV.long_name = VARS[va][3]
                    UV.units = VARS[va][4]
                    UV.timedaystat = VARS[va][2]

                    # write data to netcdf
                    lat[:]=rgrLat
                    lon[:]=rgrLon
                    UV[:]=rgrDATA
                    time[:]=time_act
                    root_grp.close()

                    # compress the netcdf file
                    pp = subprocess.Popen("nccopy -k 4 -d 1 -s "+File+' '+File_fin, shell=True)
                    pp.wait()
                    subprocess.Popen("rm  "+File, shell=True)
            pp = subprocess.Popen('rm -r '+OutDir+str(YYYYall[yy]), shell=True)
            pp.wait()


# In[35]:


FileACT


# In[ ]:




