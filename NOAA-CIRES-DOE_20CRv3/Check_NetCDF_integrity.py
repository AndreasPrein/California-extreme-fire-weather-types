#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python


# # Check_NetCDF_integrity.ipynb

# In[2]:


'''
    File name: Check_NetCDF_integrity.ipynb
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 23.07.2020
    Date last modified: 23.07.2020

    ##############################################################
    Purpos:

    1) loops over all NetCDFs in a dierctory and finds those that are not readable


'''


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
import time


# In[2]:


# ________________________________________________________________________
# ________________________________________________________________________
#             USER INPUT SECTION

DataDir='/glade/campaign/mmm/c3we/prein/NOAA-20C/daymean/'

VARS=[['PRES','PRES','mean','surface pressure','hPa'],
      ['PRMSL','PSL','mean','sea level pressure','hPa'],
      ['PWAT','PW','mean','precipitable water','mm'],
      ['SPFH2m','Q2','mean','2m specific humidity', 'g kg-1'],
      ['SPFH500','Q500','mean','500 hPa specific humidity','g kg-1'],
      ['TMP2m','T2','mean','2m mean air temperature', 'K'],
      ['TMP2m','T2min','min','2m minimum air temperature', 'K'],
      ['TMP2m','T2max','max','2m maximum air temperature', 'K'],
      ['TMP500','T500','mean','500 hPa mean air temperature', 'K'],
      ['UGRD10m','U10','mean','10m zonal wind speed', 'm s-1'],
      ['UGRD200','U200','mean','200hPa zonal wind speed', 'm s-1'],
      ['UGRD500','U500','mean','500hPa zonal wind speed', 'm s-1'],
      ['VGRD10m','V10','mean','10m meridional wind speed', 'm s-1'],
      ['VGRD200','V200','mean','200hPa meridional wind speed', 'm s-1'],
      ['VGRD500','V500','mean','500hPa meridional wind speed', 'm s-1'],
      ['VGRD850','V850','mean','850hPa meridional wind speed', 'm s-1']]

var = str(sys.argv[1])
SkipYear = int(sys.argv[2])-1900

# In[ ]:


# Loop over netcdf files in directory and check if we can read on colum of data
NCfiles = np.sort(glob.glob(DataDir+var+'/'+"*.nc"))
iStart = int(SkipYear*80)
for fi in range(len(NCfiles))[iStart:]:
    try:
        print('    '+NCfiles[fi])
        ncid=Dataset(NCfiles[fi], mode='r')
        TEST=np.squeeze(ncid.variables[var][:,100,200])
        ncid.close()
    except:
        print('    Delet corrupt file - '+NCfiles[fi])
        pp = subprocess.Popen("rm "+NCfiles[fi], shell=True); pp.wait()
    time.sleep(.3)


# In[ ]:





# In[ ]:




