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

# CA FIRE REGIONS
sPlotDir='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/plots/WT-Centroids/'# +str(iNrOfExtremes)+'_Events/'
sDataDir='/glade/campaign/mmm/c3we/prein/Papers/2019_Janice-CA-Fire-WTs/data/preselect/'
DW_Regions=['Bay_Area','Central_Coast','LA','Modoc','Northeast','San_Diego','Siearas_East','Siearas_West']
# iRegion=1
# Region=DW_Regions[iRegion]
sSubregionPR='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/Shapefiles/' # +DW_Regions[iRegion]

# rgsWTvars=['CAPE','CIN','LCL','LFC','MF2AVG','MF500','MF850','MR2AVG','MR500','MR850','PWAVG','RH2AVG','RH500','RH850','SLPAVG','T2AVG','T2MAX','T2MIN','T500','T850','U10AVG','U200','U500','U850','V10AVG','V200','V500','V850','WDIR10','WDIR200','WDIR500','WDIR850','WSPD10','WSPD200','WSPD500','WSPD850','Z500']

### Get the settings from /glade/campaign/mmm/c3we/prein/Papers/2019_Janice-CA-Fire-WTs/data/BestSettings

# # Bay Area
# rgsWTvars=['RH2AVG','T2MAX','Z500','T2AVG','MF2AVG','V850','T2MIN','MR850'] #,'MF850',V10AVG,RH500,MR500,RH850,SLPAVG,WSPD850,CAPE,WSPD10,T850,LCL,U10AVG,PWAVG,CIN,MF500,MR2AVG,V200,U850,T500,U200,WSPD200,LFC,U500,V500,WSPD500]

# # Central_Coast
# rgsWTvars=['T2AVG','T2MIN','RH2AVG','T2MAX','V10AVG','Z500','SLPAVG','RH500'] #,WSPD10,MR500,T500,RH850,MF2AVG,T850,PWAVG,MF850,V850,WSPD850,U850,MR850,WSPD500,CIN,LCL,MF500,U500,WSPD200,CAPE,MR2AVG,V500,LFC,U200,U10AVG,V200]

# # Modoc
# rgsWTvars=['U500','U10AVG','SLPAVG','T850','MR850','MR2AVG','WSPD10','V850'] #,RH2AVG,PWAVG,T2AVG,V10AVG,U850,T2MIN,T2MAX,LCL,Z500,LFC,WSPD850,U200,MF850,MF2AVG,RH850,MF500,CAPE,WSPD200,V200,WSPD500,V500,T500,MR500,CIN,RH500]

# Northeast
# rgsWTvars=['U200','T500','RH500','MR500','PWAVG','MF2AVG','T2AVG','CAPE'] #,T2MIN,T2MAX,RH2AVG,MF850,WSPD10,MR850,LFC,V10AVG,MF500,V850,MR2AVG,Z500,U10AVG,WSPD200,WSPD850,RH850,SLPAVG,T850,U850,V500,CIN,U500,V200,WSPD500,LCL]

# San_Diego
#rgsWTvars=["V200","PWAVG","MF2AVG","WSPD10","MF500","WDIR500","V10AVG","T500"] #,"V500","LCL"] ,"MR850","WDIR850","WSPD850","T2MIN","V850","U10AVG","MF850","WSPD200","U850","WDIR10","T850","SLPAVG","RH2AVG","Z500","RH850","WSPD500","U200","MR2AVG","CAPE","LFC","T2AVG","RH500","U500","WDIR200","MR500","CIN","T2MAX"]

# LA
# rgsWTvars=["CAPE","MR2AVG","LFC","SLPAVG","T850","U850","V850","MF2AVG"] #,"V10AVG","Z500","WSPD500","WSPD850","RH2AVG","WSPD10","RH850","MF850","U10AVG","MF500","T500","LCL","MR850","V500","U500","T2MIN","T2AVG","WSPD200","U200","T2MAX","RH500","V200","MR500","PWAVG","CIN"]

# # Siearra East
# rgsWTvars=['T2MIN','T2AVG','LFC','T2MAX','Z500','T500','RH850','WSPD10'] #,V10AVG,MF500,LCL,T850,U850,MR850,SLPAVG,U10AVG,MF850,MR2AVG,MF2AVG,RH500,V500,V200,WSPD200,WSPD500,WSPD850,CAPE,MR500,V850,U500,PWAVG,U200,CIN,RH2AVG]

# Sierra West
rgsWTvars=['RH2AVG','MR2AVG','V500','LFC','LCL','Z500','V200','WSPD200'] #,MR850,MF2AVG,WSPD500,RH850,CAPE,PWAVG,U200,T2MAX,V10AVG,T2MIN,MF500,V850,SLPAVG,U850,T2AVG,MF850,WSPD10,T500,U500,MR500,WSPD850,T850,CIN,RH500,U10AVG

VarsFullName   =rgsWTvars
rgsWTfolders=['/glade/campaign/mmm/c3we/prein/ERA5/' for va in range(len(rgsWTvars))]

rgrNrOfExtremes=[3,4,6,15] #[3,4,6,15,30]

WT_Domains=['M'] #['S','M','L'] # ['S','M','L','XXL'] 
DomDegreeAdd=[5] #[2,5,10]   # [2,5,10,20] 

Annual_Cycle=['1'] # '1' means that the annual cycle gets removed before clustering; '0' nothing is done

SpatialSmoothing=[0.5] #[0,0.5,1] #[0,0.5,1]

Metrics=['PSS','MRD','MRR','APR','PEX']

Dimensions=['Variables','Extreme Nr.','Domain Size','Annual Cycle','Smoothing','Split Sample','Metrics']

