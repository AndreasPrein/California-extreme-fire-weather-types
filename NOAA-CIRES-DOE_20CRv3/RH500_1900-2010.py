#!/usr/bin/env python


# # RH500_1900-2010.ipynb

# In[1]:


'''
    File name: RH500_1900-2010.ipynb
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 23.07.2020
    Date last modified: 23.07.2020

    ##############################################################
    Purpos:

    1) Reads in T2 and DT2 variables from the ERA-20C

    2) Calculate RH

    3) store the data as NetCDF


'''


# In[2]:


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
import wrf


# In[3]:


# ________________________________________________________________________
# ________________________________________________________________________
#             USER INPUT SECTION

DataDir='/glade/campaign/mmm/c3we/prein/NCEP-20C/daymean/'
SaveDir='/glade/campaign/mmm/c3we/prein/NCEP-20C/daymean/RH500/'
if not os.path.exists(SaveDir):
    os.makedirs(SaveDir)
Vars=['T500','Q500']
FinVar='RH500'
ConstantVars='/gpfs/fs1/collections/rda/data/ds131.3/anl/anl_mean_2015_VGRD_pres.nc'

rgdTime = pd.date_range(datetime.datetime(1900, 1, 1,0), end=datetime.datetime(2015, 12, 31,23), freq='d')
YYYYall = np.unique(rgdTime.year)


# In[4]:


# ________________________________________________________________________
# ________________________________________________________________________
#              READ IN THE COORDINATES
ncid=Dataset(ConstantVars, mode='r')
rgrLat=np.squeeze(ncid.variables['latitude'][:])
rgrLon=np.squeeze(ncid.variables['longitude'][:])
ncid.close()


# In[ ]:


for yy in range(len(YYYYall)):
    print('    process '+str(YYYYall[yy]))
    for me in range(80):
        T_file=DataDir+Vars[0]+'/'+Vars[0]+'.'+str(YYYYall[yy])+'_mem'+str(me+1).zfill(3)+'.nc'
        Q_file=DataDir+Vars[1]+'/'+Vars[1]+'.'+str(YYYYall[yy])+'_mem'+str(me+1).zfill(3)+'.nc'
        T_Filename=T_file.split('/')[-1]
        File_fin=SaveDir+'/'+T_Filename.replace(Vars[0], FinVar)
        File=SaveDir+T_Filename.replace(Vars[0], FinVar)+'_copy'
        if os.path.exists(File_fin) == 0:
            if not os.path.exists(SaveDir):
                os.makedirs(SaveDir)
#             print( '    '+File_fin)
            ncid=Dataset(T_file, mode='r')
            rgrT=np.squeeze(ncid.variables[Vars[0]][:])
            time_act=np.squeeze(ncid.variables['time'][:])
            ncid.close()
            ncid=Dataset(Q_file, mode='r')
            rgrQ=np.squeeze(ncid.variables[Vars[1]][:])
            ncid.close()
            
            rgrP = np.copy(rgrQ); rgrP[:] = 50000.
            RH_vals = wrf.rh(rgrQ, rgrP , rgrT, meta=False)

            # ________________________________________________________________________
            # write the netcdf
#             print( '        ----------------------')
#             print( '        Save data to '+File_fin)
            root_grp = Dataset(File, 'w', format='NETCDF4')
            # dimensions
            root_grp.createDimension('time', None)
            root_grp.createDimension('rlon', rgrLon.shape[0])
            root_grp.createDimension('rlat', rgrLat.shape[0])
            # variables
            lat = root_grp.createVariable('lat', 'f4', ('rlat',))
            lon = root_grp.createVariable('lon', 'f4', ('rlon',))
            time = root_grp.createVariable('time', 'f8', ('time',))
            RH = root_grp.createVariable(FinVar, 'f4', ('time','rlat','rlon',),fill_value=-99999)

            time.calendar = "proleptic_gregorian"
            time.units = "hours since 1900-1-1 00:00:00"
            time.standard_name = "time"
            time.long_name = "time"

            lon.standard_name = "longitude"
            lon.long_name = "longitude"
            lon.units = "degrees_east"

            lat.standard_name = "latitude"
            lat.long_name = "latitude"
            lat.units = "degrees_north"

            RH.standard_name = "500 hPa relative humidity"
            RH.long_name = "500 hPa relative humidity"
            RH.units = "%"

            # write data to netcdf
            lat[:]=rgrLat
            lon[:]=rgrLon
            RH[:]=RH_vals
            time[:]=time_act
            root_grp.close()

            # compress the netcdf file
            pp = subprocess.Popen("nccopy -k 4 -d 1 -s "+File+' '+File_fin, shell=True)
            pp.wait()
            subprocess.Popen("rm  "+File, shell=True)

