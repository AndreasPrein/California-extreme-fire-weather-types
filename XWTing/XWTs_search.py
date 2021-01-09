#!/usr/bin/env python
'''File name: Denver-Water_XWT.py
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 16.04.2018
    Date last modified: 16.04.2018

    ############################################################## 
    Purpos:
    Contains the setup for extreme weather typing (XWT) for
    Denver Water watersheds

'''

from pdb import set_trace as stop
import numpy as np
import os
import pandas as pd
import datetime

# ###################################################

FireObs='Parks' # ['MODIS','Parks']
if FireObs == 'MODIS':
    dStartDayPR=datetime.datetime(2001, 01, 01,0) # (2000, 11, 01,0)
    dStopDayPR=datetime.datetime(2015, 12, 31,23)
    rgdTime = pd.date_range(dStartDayPR, end=dStopDayPR, freq='d')
elif FireObs == 'Parks':
    dStartDayPR=datetime.datetime(2002, 01, 01,0) # (2000, 11, 01,0)
    dStopDayPR=datetime.datetime(2018, 12, 31,23)
    rgdTime = pd.date_range(dStartDayPR, end=dStopDayPR, freq='d')
rgdTime = pd.date_range(dStartDayPR, end=dStopDayPR, freq='d')
iMonths=[1,2,3,4,5,6,7,8,9,10,11,12] # [1,2,3,10,11,12] # [4,5,6,7,8,9]

# CA FIRE REGIONS
sPlotDir='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/plots/WT-Centroids/'# +str(iNrOfExtremes)+'_Events/'
sDataDir='/glade/campaign/mmm/c3we/prein/Papers/2019_Janice-CA-Fire-WTs/data/'
DW_Regions=['Bay_Area','Central_Coast','LA','Modoc','Northeast','San_Diego','Siearas_East','Siearas_West']
# iRegion=1
# Region=DW_Regions[iRegion]
sSubregionPR='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/Shapefiles/' # +DW_Regions[iRegion]

rgsWTvars= ['var151','u',   'v',     'UV', 'tcw',  'FLX',    'q',    'q'   , 'z','cape', 'UV']
VarsFullName=['PSL','U850','V850','UV850',  'PW',  'FLX850', 'Q850', 'Q500', 'ZG500','CAPE', 'UV200']
rgsWTfolders=['/glade/scratch/prein/ERA-Interim/PSL/fin_PSL-sfc_ERA-Interim_12-0_',\
              '/glade/scratch/prein/ERA-Interim/UV850/fin_UV850-sfc_ERA-Interim__',\
              '/glade/scratch/prein/ERA-Interim/UV850/fin_UV850-sfc_ERA-Interim__',\
              '/glade/scratch/prein/ERA-Interim/UV850/fin_UV850-sfc_ERA-Interim__',\
              '/glade/scratch/prein/ERA-Interim/TCW/fin_TCW-sfc_ERA-Interim_12-0_',\
              '/glade/scratch/prein/ERA-Interim/UV850xTCW/fin_FLX-pl_ERA-Interim_',\
              '/glade/scratch/prein/ERA-Interim/Q850/Q850_daymean_',\
              '/glade/scratch/prein/ERA-Interim/Q500/Q500_daymean_',\
              '/glade/scratch/prein/ERA-Interim/Z500/Z500_daymean_',\
              '/glade/scratch/prein/ERA-Interim/CAPE_ECMWF/fin_CAPE-ECMWF-sfc_ERA-Interim_12-9_',\
              '/glade/scratch/prein/ERA-Interim/UV200/UV200_daymean_']



rgrNrOfExtremes=[3,4,6,15] #[3,4,6,15,30]

WT_Domains=['L'] #['S','M','L'] # ['S','M','L','XXL'] 
DomDegreeAdd=[10] #[2,5,10]   # [2,5,10,20] 

Annual_Cycle=['1'] # '1' means that the annual cycle gets removed before clustering; '0' nothing is done

SpatialSmoothing=[0.5] #[0,0.5,1] #[0,0.5,1]

Metrics=['PSS','MRD','MRR','APR','PEX']

Dimensions=['Variables','Extreme Nr.','Domain Size','Annual Cycle','Smoothing','Split Sample','Metrics']

