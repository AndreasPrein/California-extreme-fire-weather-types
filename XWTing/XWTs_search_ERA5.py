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

FireObs='MODIS' # ['MODIS','Parks']
if FireObs == 'MODIS':
    dStartDayPR=datetime.datetime(2001, 1, 1,0) # (2001, 01, 01,0)
    dStopDayPR=datetime.datetime(2019, 11, 30,23)
    rgdTime = pd.date_range(dStartDayPR, end=dStopDayPR, freq='d')
elif FireObs == 'Parks':
    dStartDayPR=datetime.datetime(2002, 1, 1,0) # (2000, 11, 01,0)
    dStopDayPR=datetime.datetime(2018, 12, 31,23)
    rgdTime = pd.date_range(dStartDayPR, end=dStopDayPR, freq='d')
rgdTime = pd.date_range(dStartDayPR, end=dStopDayPR, freq='d')
iMonths=[1,2,3,4,5,6,7,8,9,10,11,12] # [1,2,3,10,11,12] # [4,5,6,7,8,9]

# ---------
# Setup clustering algorithm
ClusterMeth='hdbscan'  # current options are ['HandK','hdbscan']
ClusterBreakup = 0     # breakes up clusters that are unproportionally large (only for hdbscan)
RelAnnom=1           # 1 - calculates daily relative anomalies
NormalizeData='D'    # normalize variables | options are  - 'C' - climatology
                                                        # - 'D' - daily (default)
                                                        # - 'N' - none
MinDistDD=7          # minimum nr of days between XWT events
RemoveAnnualCycl=1   # remove annual cycle in varaiables with 21 day moving average filter
# ---------


# CA FIRE REGIONS
sPlotDir='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/plots/WT-Centroids/'# +str(iNrOfExtremes)+'_Events/'
sDataDir='/glade/campaign/mmm/c3we/prein/Papers/2019_Janice-CA-Fire-WTs/data/AUC-APR/'
DW_Regions=['Bay_Area','Central_Coast','LA','Modoc','Northeast','San_Diego','Siearas_East','Siearas_West']
# iRegion=1
# Region=DW_Regions[iRegion]
sSubregionPR='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/Shapefiles/' # +DW_Regions[iRegion]

# rgsWTvars=['MF2AVG','MF500','MF850','MR2AVG','MR500','MR850','PWAVG','RH2AVG','RH500','RH850','SLPAVG','T2AVG','T2MAX','T2MIN','T500','T850','U10AVG','U200','U500','U850','V10AVG','V200','V500','V850','WSPD10','WSPD200','WSPD500','WSPD850','Z500']

rgsWTvars=['CAPE','CIN','LCL','LFC','MF2AVG','MF500','MF850','MR2AVG','MR500','MR850','PWAVG','RH2AVG','RH500','RH850','SLPAVG','T2AVG','T2MAX','T2MIN','T500','T850','U10AVG','U200','U500','U850','V10AVG','V200','V500','V850','WSPD10','WSPD200','WSPD500','WSPD850','Z500']



VarsFullName   = rgsWTvars
rgsWTfolders=['/glade/campaign/mmm/c3we/prein/ERA5/' for va in range(len(rgsWTvars))]

rgrNrOfExtremes=[4,6,10,15] #,6,10,14] #[3,4,6,15] #[3,4,6,15,30]

WT_Domains=['S','M'] #['S','M','L'] # ['S','M','L','XXL']
DomDegreeAdd=[2,5] #[2,5,10]   # [2,5,10,20] 

Annual_Cycle=['1'] # '1' means that the annual cycle gets removed before clustering; '0' nothing is done

SpatialSmoothing=[0,0.5] #[0,0.5,1] #[0,0.5,1]

Metrics=['PSS','MRD','MRR','APR','PEX','AUC']

Dimensions=['Variables','Extreme Nr.','Domain Size','Annual Cycle','Smoothing','Split Sample','Metrics']

