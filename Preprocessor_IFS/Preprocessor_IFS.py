#!/usr/bin/env python


# # Preprocessor_IFS.ipynb

# In[1]:


'''
    File name: Preprocessor_IFS.ipynb
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 22.09.2020
    Date last modified: 22.09.2020

    ##############################################################
    Purpos:
    
    Reads in the IFS data from RDA and saves annual files over focus
    region that can be used in the XWT program


'''


# In[20]:


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
from thermodynamics_p3 import TandDTtoRH, VaporPressure, MixRatio, relhum
import wrf


def WriteNetCDF(File_fin,
               File,
               LON,
               LAT,
               rgdTimeMM,
               TIME,
               VarName,
               VarLongName,
               VarUnit,
               DATA):


    # ________________________________________________________________________
    # write the netcdf
    print( '        ----------------------')
    print( '        Save data to '+File_fin)
    root_grp = Dataset(File, 'w', format='NETCDF4')
    # dimensions
    root_grp.createDimension('time', None)
    root_grp.createDimension('longitude', LON.shape[0])
    root_grp.createDimension('latitude', LAT.shape[0])
    # variables
    lat = root_grp.createVariable('latitude', 'f4', ('latitude',))
    lon = root_grp.createVariable('longitude', 'f4', ('longitude',))
    time = root_grp.createVariable('time', 'f8', ('time',))
    MF = root_grp.createVariable(VarName, 'f4', ('time','latitude','longitude',),fill_value=-99999)

    time.calendar = "proleptic_gregorian"
    time.units = "days since "+str(rgdTimeMM[0].year)+"-"+str(rgdTimeMM[0].month)+"-1 12:00:00"
    time.standard_name = "time"
    time.long_name = "time"

    lon.standard_name = "longitude"
    lon.long_name = "longitude"
    lon.units = "degrees_east"
    lon.short_name = "lon"

    lat.standard_name = "latitude"
    lat.long_name = "latitude"
    lat.units = "degrees_north"
    lat.qualifier = "Gaussian"
    lat.short_name = "lat"

    MF.standard_name = VarLongName
    MF.long_name = VarLongName
    MF.units = VarUnit

    # write data to netcdf
    lat[:]=LAT
    lon[:]=LON
    MF[:]=DATA
    time[:]=TIME
    root_grp.close()

    # compress the netcdf file
    pp = subprocess.Popen("nccopy -k 4 -d 1 -s "+File+' '+File_fin, shell=True)
    pp.wait()
    subprocess.Popen("rm  "+File, shell=True)


# In[21]:


# ________________________________________________________________________
# ________________________________________________________________________
#             USER INPUT SECTION

YYYY = int(sys.argv[1]) #2020

DataDirPL = '/gpfs/fs1/collections/rda/data/ds113.1/ec.oper.an.pl/'   # ec.oper.an.pl.128_060_pv.regn1280sc.2020080100.nc
DataDirSL = '/gpfs/fs1/collections/rda/data/ds113.1/ec.oper.fc.sfc/'   # 202008/ec.oper.fc.sfc.128_031_ci.regn1280sc.20200801.nc
SaveDir='/glade/campaign/mmm/c3we/prein/IFS/'
VarsSL =   [['T2AVG','167_2t','VAR_2T','2 metre temperature',"K"],
            ['T2MIN','202_mn2t','MN2T','Minimum temperature at 2 metres','K'],
            ['T2MAX','201_mx2t','MX2T','Maximum temperature at 2 metres','K'],
            ['PWAVG','136_tcw','TCW','Precipitable Water',"kg m**-2"],
            ['SLPAVG','151_msl','MSL',"Mean sea level pressure",'Pa'],
            ['SP','134_sp','SP','Surface pressure','Pa'],
            ['DT2','168_2d','VAR_2D',"2 metre dewpoint temperature",'K'],
            ['U10','165_10u','VAR_10U',"10 metre U wind component","m s**-1"],
            ['V10','166_10v','VAR_10V',"10 metre V wind component","m s**-1"]]

VarsPL = [['T','130_t','T',"Temperature",'K'],
          ['V','132_v','V',"V component of wind","m s**-1"],
          ['U','131_u','U',"U component of wind","m s**-1"],
          ['Q','133_q','Q',"Specific humidity","kg kg**-1"],
          ['Z','129_z',"Z","Geopotential","m**2 s**-2"]]

Plevs=['850','500','200']
HH = ['00','06','12','18']
FocusReg = [50,360-96,25,360-140]

rgdTimeDD = pd.date_range(datetime.datetime(2016, 1, 1,0), end=datetime.datetime(2020, 11, 30,23), freq='d')
rgiTimeDD = np.array(range(len(rgdTimeDD)))
rgdTimeMM = pd.date_range(datetime.datetime(2016, 1, 1,0), end=datetime.datetime(2020, 11, 30,23), freq='m')
YYYYall = np.unique(rgdTimeDD.year)

MONrange = np.array(range(len(rgdTimeMM)))[rgdTimeMM.month == YYYY]


# In[22]:


# ________________________________________________________________________
# ________________________________________________________________________
#              READ IN THE COORDINATES
ncid=Dataset('/gpfs/fs1/collections/rda/data/ds113.1/ec.oper.an.pl/202008/ec.oper.an.pl.128_133_q.regn1280sc.2020082418.nc', mode='r')
rgrLat=np.squeeze(ncid.variables['latitude'][:])
rgrLon=np.squeeze(ncid.variables['longitude'][:])
Plev = np.squeeze(ncid.variables['level'][:])
ncid.close()

# Find region of interest
NN = np.where(rgrLat >= FocusReg[0])[0][-1]
EE = np.where(rgrLon >= FocusReg[1])[0][0]
SS = np.where(rgrLat <= FocusReg[2])[0][0]
WW = np.where(rgrLon <= FocusReg[3])[0][-1]

LAT = rgrLat[NN:SS]
LON = rgrLon[WW:EE]

iPL = np.array([np.where(Plev == int(Plevs[pl]))[0][0] for pl in range(len(Plevs))])


# ### Work on the surface level files

# In[23]:

for mm in MONrange:
    for sl in range(len(VarsSL)):
        print('    process '+str(VarsSL[sl][0]))
        SaveDirFIN = SaveDir+VarsSL[sl][0]
        sYYYYMM = str(rgdTimeMM[mm].year)+str(rgdTimeMM[mm].month).zfill(2)
        File_fin=SaveDirFIN+'/'+VarsSL[sl][0]+'_'+sYYYYMM+'_daymean.nc'
        File=File_fin+'_copy'
        if os.path.exists(File_fin) == 0:
            if not os.path.exists(SaveDirFIN):
                os.makedirs(SaveDirFIN)
            print('        '+str(rgdTimeMM[mm]))
            DATA = np.zeros((rgdTimeMM[mm].day,len(LAT),len(LON)))
            TIME = rgiTimeDD[(rgdTimeDD.month == rgdTimeMM[mm].month) & (rgdTimeDD.year == rgdTimeMM[mm].year)]
            for dd in range(rgdTimeMM[mm].day):
                sYYYYMMDD = str(rgdTimeMM[mm].year)+str(rgdTimeMM[mm].month).zfill(2)+str(dd+1).zfill(2)
                DATAfile = DataDirSL+sYYYYMM+'/ec.oper.fc.sfc.128_'+VarsSL[sl][1]+'.regn1280sc.'+sYYYYMMDD+'.nc'
                ncid=Dataset(DATAfile, mode='r')
                DATA[dd,:,:] = np.mean(np.squeeze(ncid.variables[VarsSL[sl][2]][:,:,NN:SS,WW:EE]), axis=(0,1))
                ncid.close()

            WriteNetCDF(File_fin,
                   File,
                   LON,
                   LAT,
                   rgdTimeMM,
                   TIME,
                   VarsSL[sl][0],
                   VarsSL[sl][3],
                   VarsSL[sl][4],
                   DATA)


# ### Work on pressure level data

# In[4]:


for mm in MONrange:
    for pl in range(len(Plevs)):
        for sl in range(len(VarsPL)):
            print('    process '+str(VarsPL[sl][0]))
            SaveDirFIN = SaveDir+VarsPL[sl][0]+Plevs[pl]
            sYYYYMM = str(rgdTimeMM[mm].year)+str(rgdTimeMM[mm].month).zfill(2)
            File_fin=SaveDirFIN+'/'+VarsPL[sl][0]+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
            File=File_fin+'_copy'
            if os.path.exists(File_fin) == 0:
                if not os.path.exists(SaveDirFIN):
                    os.makedirs(SaveDirFIN)
                print('        '+str(rgdTimeMM[mm]))
                DATA = np.zeros((rgdTimeMM[mm].day,len(HH),len(LAT),len(LON)))
                TIME = rgiTimeDD[(rgdTimeDD.month == rgdTimeMM[mm].month) & (rgdTimeDD.year == rgdTimeMM[mm].year)]
                for dd in range(rgdTimeMM[mm].day):
                    sYYYYMMDD = str(rgdTimeMM[mm].year)+str(rgdTimeMM[mm].month).zfill(2)+str(dd+1).zfill(2)
                    print('            '+sYYYYMMDD)
                    for hh in range(len(HH)):
                        if ((VarsPL[sl][0] != 'U') & (VarsPL[sl][0] != 'V')):
                            DATAfile = DataDirPL+sYYYYMM+'/ec.oper.an.pl.128_'+VarsPL[sl][1]+'.regn1280sc.'+sYYYYMMDD+HH[hh]+'.nc'
                        else:
                            DATAfile = DataDirPL+sYYYYMM+'/ec.oper.an.pl.128_'+VarsPL[sl][1]+'.regn1280uv.'+sYYYYMMDD+HH[hh]+'.nc'
                        ncid=Dataset(DATAfile, mode='r')
                        DATA[dd,:,:] = np.squeeze(ncid.variables[VarsPL[sl][2]][:,iPL[pl],NN:SS,WW:EE])
                        ncid.close()
                        
                        
                WriteNetCDF(File_fin,
                   File,
                   LON,
                   LAT,
                   rgdTimeMM,
                   TIME,
                   VarsPL[sl][0]+Plevs[pl],
                   VarsPL[sl][3]+' at '+Plevs[pl]+' hPa',
                   VarsPL[sl][4],
                   np.mean(DATA, axis=1))
                        
#                 # ________________________________________________________________________
#                 # write the netcdf
#                 print( '        ----------------------')
#                 print( '        Save data to '+File_fin)
#                 root_grp = Dataset(File, 'w', format='NETCDF4')
#                 # dimensions
#                 root_grp.createDimension('time', None)
#                 root_grp.createDimension('longitude', LON.shape[0])
#                 root_grp.createDimension('latitude', LAT.shape[0])
#                 # variables
#                 lat = root_grp.createVariable('latitude', 'f4', ('latitude',))
#                 lon = root_grp.createVariable('longitude', 'f4', ('longitude',))
#                 time = root_grp.createVariable('time', 'f8', ('time',))
#                 MF = root_grp.createVariable(VarsPL[sl][0]+Plevs[pl], 'f4', ('time','latitude','longitude',),fill_value=-99999)

#                 time.calendar = "proleptic_gregorian"
#                 time.units = "days since "+str(rgdTimeMM[0].year)+"-"+str(rgdTimeMM[0].month)+"-1 12:00:00"
#                 time.standard_name = "time"
#                 time.long_name = "time"

#                 lon.standard_name = "longitude"
#                 lon.long_name = "longitude"
#                 lon.units = "degrees_east"
#                 lon.short_name = "lon"

#                 lat.standard_name = "latitude"
#                 lat.long_name = "latitude"
#                 lat.units = "degrees_north"
#                 lat.qualifier = "Gaussian"
#                 lat.short_name = "lat"

#                 MF.standard_name = VarsPL[sl][3]+' at '+Plevs[pl]+' hPa'
#                 MF.long_name = VarsPL[sl][3]+' at '+Plevs[pl]+' hPa'
#                 MF.units = VarsPL[sl][4]

#                 # write data to netcdf
#                 lat[:]=LAT
#                 lon[:]=LON
#                 MF[:]= np.mean(DATA, axis=1)
#                 time[:]=TIME
#                 root_grp.close()

#                 # compress the netcdf file
#                 pp = subprocess.Popen("nccopy -k 4 -d 1 -s "+File+' '+File_fin, shell=True)
#                 pp.wait()
#                 subprocess.Popen("rm  "+File, shell=True)


# ### Calculate combined variables

# In[11]:


for mm in MONrange:
    sYYYYMM = str(rgdTimeMM[mm].year)+str(rgdTimeMM[mm].month).zfill(2)
    TIME = rgiTimeDD[(rgdTimeDD.month == rgdTimeMM[mm].month) & (rgdTimeDD.year == rgdTimeMM[mm].year)]
    
    # Start with 10m wind fields
    SaveDirFIN = SaveDir+'UV10'
    File_fin=SaveDirFIN+'/UV10_'+sYYYYMM+'_daymean.nc'
    File=File_fin+'_copy'
    if os.path.exists(File_fin) == 0:
        if not os.path.exists(SaveDirFIN):
            os.makedirs(SaveDirFIN)
    
        Uf = SaveDir+'U10/U10'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Uf, mode='r')
        U10=np.squeeze(ncid.variables['U10'][:])
        ncid.close()
        Vf = SaveDir+'V10/V10'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Vf, mode='r')
        V10=np.squeeze(ncid.variables['V10'][:])
        ncid.close()
        UV10 = (U10**2+V10**2)**0.5

        WriteNetCDF(File_fin,
                   File,
                   LON,
                   LAT,
                   rgdTimeMM,
                   TIME,
                   'UV10',
                   'wind speed 10 m above surface',
                   'm s**-2',
                   UV10)
    # Windspeed on pressure levels
    for pl in range(len(Plevs)):
        SaveDirFIN = SaveDir+'UV'+Plevs[pl]
        File_fin=SaveDirFIN+'/'+'UV'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
        File=File_fin+'_copy'
        if os.path.exists(File_fin) == 0:
            if not os.path.exists(SaveDirFIN):
                os.makedirs(SaveDirFIN)
            Uf = SaveDir+'U'+Plevs[pl]+'/U'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
            ncid=Dataset(Uf, mode='r')
            UU=np.squeeze(ncid.variables['U'+Plevs[pl]][:])
            ncid.close()
            Vf = SaveDir+'V'+Plevs[pl]+'/V'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
            ncid=Dataset(Vf, mode='r')
            VV=np.squeeze(ncid.variables['V'+Plevs[pl]][:])
            ncid.close()
            UV = (UU**2+VV**2)**0.5

            WriteNetCDF(File_fin,
                       File,
                       LON,
                       LAT,
                       rgdTimeMM,
                       TIME,
                       'UV'+Plevs[pl],
                       'wind speed at '+Plevs[pl],
                       'm s**-2',
                       UV)
    
    # Calculate Q at 2 m
    SaveDirFIN = SaveDir+'Q2'
    File_fin=SaveDirFIN+'/Q2_'+sYYYYMM+'_daymean.nc'
    File=File_fin+'_copy'
    if os.path.exists(File_fin) == 0:
        if not os.path.exists(SaveDirFIN):
            os.makedirs(SaveDirFIN)
    
        Uf = SaveDir+'T2AVG/T2AVG'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Uf, mode='r')
        rgrT=np.squeeze(ncid.variables['T2AVG'][:])
        ncid.close()
        Vf = SaveDir+'DT2/DT2'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Vf, mode='r')
        rgrDT=np.squeeze(ncid.variables['DT2'][:])
        ncid.close()
        Vf = SaveDir+'SP/SP'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Vf, mode='r')
        rgrSP=np.squeeze(ncid.variables['SP'][:])
        ncid.close()
        
        RH_vals=TandDTtoRH(rgrT,rgrDT)
        Ws = VaporPressure(rgrDT-273.15)
        W = (RH_vals/100.) * Ws
        Q2 = MixRatio(W,rgrSP)

        WriteNetCDF(File_fin,
                   File,
                   LON,
                   LAT,
                   rgdTimeMM,
                   TIME,
                   'Q2',
                   'mixing ratio 2m above surface',
                   'km kg**-1',
                   Q2)
        
    # Calculate 2m moisture flux
    SaveDirFIN = SaveDir+'MF2AVG'
    File_fin=SaveDirFIN+'/MF2AVG_'+sYYYYMM+'_daymean.nc'
    File=File_fin+'_copy'
    if os.path.exists(File_fin) == 0:
        if not os.path.exists(SaveDirFIN):
            os.makedirs(SaveDirFIN)
    
        Uf = SaveDir+'Q2/Q2'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Uf, mode='r')
        QQ=np.squeeze(ncid.variables['Q2'][:])
        ncid.close()
        Vf = SaveDir+'UV10/UV10'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Vf, mode='r')
        UV=np.squeeze(ncid.variables['UV10'][:])
        ncid.close()
        MF = QQ*UV

        WriteNetCDF(File_fin,
                   File,
                   LON,
                   LAT,
                   rgdTimeMM,
                   TIME,
                   'MF2AVG',
                   'moisture flux 2m above surface',
                   'km kg**-1 m s**-1',
                   MF)
        
    # Calcualte pressure level moisture flux
    for pl in range(len(Plevs)):
        SaveDirFIN = SaveDir+'FLX'+Plevs[pl]
        File_fin=SaveDirFIN+'/'+'FLX'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
        File=File_fin+'_copy'
        if os.path.exists(File_fin) == 0:
            if not os.path.exists(SaveDirFIN):
                os.makedirs(SaveDirFIN)
            Uf = SaveDir+'UV'+Plevs[pl]+'/UV'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
            ncid=Dataset(Uf, mode='r')
            UV=np.squeeze(ncid.variables['UV'+Plevs[pl]][:])
            ncid.close()
            Vf = SaveDir+'Q'+Plevs[pl]+'/Q'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
            ncid=Dataset(Vf, mode='r')
            QQ=np.squeeze(ncid.variables['Q'+Plevs[pl]][:])
            ncid.close()
            MF = QQ*UV

            WriteNetCDF(File_fin,
                       File,
                       LON,
                       LAT,
                       rgdTimeMM,
                       TIME,
                       'FLX'+Plevs[pl],
                       'moisture flux at '+Plevs[pl],
                       'km kg**-1 m s**-1',
                       MF)
            
    # Calcualte RH2m
    SaveDirFIN = SaveDir+'RH2'
    File_fin=SaveDirFIN+'/RH2_'+sYYYYMM+'_daymean.nc'
    File=File_fin+'_copy'
    if os.path.exists(File_fin) == 0:
        if not os.path.exists(SaveDirFIN):
            os.makedirs(SaveDirFIN)
    
        Uf = SaveDir+'T2AVG/T2AVG'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Uf, mode='r')
        rgrT=np.squeeze(ncid.variables['T2AVG'][:])
        ncid.close()
        Vf = SaveDir+'DT2/DT2'+'_'+sYYYYMM+'_daymean.nc'
        ncid=Dataset(Vf, mode='r')
        rgrDT=np.squeeze(ncid.variables['DT2'][:])
        ncid.close()
        
        RH_vals=TandDTtoRH(rgrT,rgrDT)

        WriteNetCDF(File_fin,
                   File,
                   LON,
                   LAT,
                   rgdTimeMM,
                   TIME,
                   'RH2',
                   'relative humidity 2m above ground',
                   '%',
                   RH_vals)
        
        
    # Calcualte RH on pressure level
    for pl in range(len(Plevs)):
        SaveDirFIN = SaveDir+'RH'+Plevs[pl]
        File_fin=SaveDirFIN+'/'+'RH'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
        File=File_fin+'_copy'
        if os.path.exists(File_fin) == 0:
            if not os.path.exists(SaveDirFIN):
                os.makedirs(SaveDirFIN)
            Uf = SaveDir+'T'+Plevs[pl]+'/T'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
            ncid=Dataset(Uf, mode='r')
            TT=np.squeeze(ncid.variables['T'+Plevs[pl]][:])
            ncid.close()
            Vf = SaveDir+'Q'+Plevs[pl]+'/Q'+Plevs[pl]+'_'+sYYYYMM+'_daymean.nc'
            ncid=Dataset(Vf, mode='r')
            QQ=np.squeeze(ncid.variables['Q'+Plevs[pl]][:])
            ncid.close()
            
            rgrP = np.copy(QQ); rgrP[:] = int(Plevs[pl])*100.
            RH_vals = wrf.rh(QQ, rgrP , TT, meta=False)

            WriteNetCDF(File_fin,
                       File,
                       LON,
                       LAT,
                       rgdTimeMM,
                       TIME,
                       'RH'+Plevs[pl],
                       'relative humidity '+Plevs[pl],
                       '%',
                       RH_vals)

