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

def HUC2_XWTs_apply(Season,
                    Region):

    from pdb import set_trace as stop
    import numpy as np
    import os
    import pandas as pd
    import datetime
    
    # ###################################################
    
    FireObs='MODIS' # ['MODIS','Parks']
    REANAL='ERA5' # ['ERAI','ERA5']
    if FireObs == 'MODIS':
        dStartDayPR=datetime.datetime(2001, 1, 1,0) # (2000, 11, 01,0)
        dStopDayPR=datetime.datetime(2019, 11, 30,23)
        rgdTime = pd.date_range(dStartDayPR, end=dStopDayPR, freq='d')
    elif FireObs == 'Parks':
        dStartDayPR=datetime.datetime(2002, 1, 1,0) # (2000, 11, 01,0)
        dStopDayPR=datetime.datetime(2018, 12, 31,23)
        rgdTime = pd.date_range(dStartDayPR, end=dStopDayPR, freq='d')
        
    if Season == 'AMJJAS':
        iMonths=[4,5,6,7,8,9]
    elif Season == 'ONDJFM':
        iMonths=[1,2,3,10,11,12]
    elif Season == 'Annual':
        iMonths=[1,2,3,4,5,6,7,8,9,10,11,12]
        
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

    # DENVER WATER REGIONS
    sPlotDir='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/plots/WT-Centroids/'# +str(iNrOfExtremes)+'_Events/'
    sDataDir='/glade/campaign/mmm/c3we/prein/Papers/2019_Janice-CA-Fire-WTs/data/AUC-APR/' # '/glade/work/jaye/fire/'
    DW_Regions=['Bay_Area','Central_Coast','LA','Modoc','Northeast','San_Diego','Siearas_East','Siearas_West']
    # sRegion=Regions.index(Region)
    # Region=Regions[sRegion]
    sSubregionPR='/glade/u/home/prein/projects/2019_Janice-CA-Fire-WTs/Shapefiles/' #+Regions[sRegion]
    
    Metrics=['PSS','MRD','MRR','APR','PEX']
    Dimensions=['Variables','Extreme Nr.','Domain Size','Annual Cycle','Smoothing','Split Sample','Metrics']


    if (Season == 'Annual') & (Region == 'Bay_Area'):
        VarsFullName=['RH2AVG','MF2AVG','V850'] #["RH2AVG","RH500","SLPAVG"] 
        rgrNrOfExtremes=4 #[6,10,15,30]
        WT_Domains='M'    # ['S','M','L','XXL'] 
        Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
        SpatialSmoothing=0.5 #[0,0.5,1]
#         VarsFullName=['RH2AVG','MF2AVG','V850'] #['MF850','MR500','RH2AVG'] #['RH2AVG','MF2AVG','V850']
#         rgrNrOfExtremes=4 #[6,10,15,30]
#         WT_Domains='M'    # ['S','M','L','XXL'] 
#         Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
#         SpatialSmoothing=0.5 #[0,0.5,1]
    if (Season == 'Annual') & (Region == 'Central_Coast'):
        VarsFullName=['RH2AVG','SLPAVG','WSPD10'] #['RH2AVG','T2MIN','WSPD10'] # ['T2AVG','T2MIN','RH2AVG']  #['RH500','SLPAVG','V10AVG']
        rgrNrOfExtremes=6 #[6,10,15,30]
        WT_Domains='M'    # ['S','M','L','XXL'] 
        Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
        SpatialSmoothing=0.5 #[0,0.5,1]
#         VarsFullName=['RH500','SLPAVG','V10AVG'] # ['T2AVG','T2MIN','RH2AVG']  #['RH500','SLPAVG','V10AVG']
#         rgrNrOfExtremes=4 #[6,10,15,30]
#         WT_Domains='M'    # ['S','M','L','XXL'] 
#         Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
#         SpatialSmoothing=0.5 #[0,0.5,1]
    if (Season == 'Annual') & (Region == 'LA'):
        VarsFullName=['MR2AVG','SLPAVG','WSPD10'] #['MR2AVG','RH2AVG','Z500']
        rgrNrOfExtremes=10 #8
        WT_Domains='M'    # ['S','M','L','XXL'] 
        Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
        SpatialSmoothing=0.5 #0
#         VarsFullName=['MR2AVG','SLPAVG','WSPD10'] # ['CAPE','MF500','RH2AVG'] #['V850','MF2AVG']
#         rgrNrOfExtremes=10 #10 #[6,10,15,30]
#         WT_Domains='M'    # ['S','M','L','XXL'] 
#         Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
#         SpatialSmoothing=0.5 #[0,0.5,1]
    if (Season == 'Annual') & (Region == 'Modoc'):
        VarsFullName=['T2MIN','WSPD10','Z500'] #['T2AVG','T2MAX','U10AVG']
        rgrNrOfExtremes=8 #[6,10,15,30]
        WT_Domains='M'    # ['S','M','L','XXL'] 
        Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
        SpatialSmoothing=0.5 #[0,0.5,1]
#         VarsFullName=['T2MIN','U10AVG','U500'] #['U500','T850','MR850']
#         rgrNrOfExtremes=8 #3 #[6,10,15,30]
#         WT_Domains='M'    # ['S','M','L','XXL'] 
#         Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
#         SpatialSmoothing=0.5 #[0,0.5,1]
    if (Season == 'Annual') & (Region == 'Northeast'):
        VarsFullName=['RH2AVG','RH500','U200'] #['MF850','RH500','T500']
        rgrNrOfExtremes=10 #[6,10,15,30]
        WT_Domains='M'    # ['S','M','L','XXL'] 
        Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
        SpatialSmoothing=0.5 #[0,0.5,1]
#         VarsFullName=['RH500','T2AVG','T500'] #['MF850','RH500','T500']
#         rgrNrOfExtremes=4 #[6,10,15,30]
#         WT_Domains='M'    # ['S','M','L','XXL'] 
#         Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
#         SpatialSmoothing=0.5 #[0,0.5,1]
    if (Season == 'Annual') & (Region == 'San_Diego'):
        if FireObs == 'MODIS':
            VarsFullName=['MF2AVG','PWAVG','V200'] # ['MF2AVG','MF500','SLPAVG']
            rgrNrOfExtremes=4 #6 #3 #[6,10,15,30]
            WT_Domains='M'    # ['S','M','L','XXL'] 
            Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
            SpatialSmoothing=0.5 #[0,0.5,1]
#             VarsFullName=['MF2AVG','PWAVG','V200']
#             rgrNrOfExtremes=4 #6 #3 #[6,10,15,30]
#             WT_Domains='M'    # ['S','M','L','XXL'] 
#             Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
#             SpatialSmoothing=0.5 #[0,0.5,1]
        elif FireObs == 'Parks':
            VarsFullName=['SLPAVG','UV200']
            rgrNrOfExtremes=6 #6 #3 #[6,10,15,30]
            WT_Domains='L'    # ['S','M','L','XXL'] 
            Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
            SpatialSmoothing=0.5 #[0,0.5,1]
    if (Season == 'Annual') & (Region == 'Siearas_East'):
        VarsFullName=['RH500','T2MIN','V200'] #['LFC','T2MIN','U10AVG']
        rgrNrOfExtremes=10 #[6,10,15,30]
        WT_Domains='M'    # ['S','M','L','XXL'] 
        Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
        SpatialSmoothing=0 #[0,0.5,1]
#         VarsFullName=['T2MIN','T500','WSPD500'] #['LFC','T2MIN','U10AVG']
#         rgrNrOfExtremes=15 #[6,10,15,30]
#         WT_Domains='M'    # ['S','M','L','XXL'] 
#         Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
#         SpatialSmoothing=0.5 #[0,0.5,1]
    if (Season == 'Annual') & (Region == 'Siearas_West'):
        VarsFullName=['RH2AVG','V500','WSPD200'] # ['PWAVG','RH2AVG','Z500'] 
        rgrNrOfExtremes=8 #[6,10,15,30]
        WT_Domains='M'    # ['S','M','L','XXL'] 
        Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
        SpatialSmoothing=0.5 #[0,0.5,1]
#         VarsFullName=['RH2AVG','V500','WSPD200']
#         rgrNrOfExtremes=8 #[6,10,15,30]
#         WT_Domains='M'    # ['S','M','L','XXL'] 
#         Annual_Cycle='1' # '1' means that the annual cycle gets removed before clustering; '0' nothing is done
#         SpatialSmoothing=0.5 #[0,0.5,1]


    # ---------------
    # Full list of available variables
    VarsFullNameAll=['CAPE', 'CIN','LCL','LFC','MF2AVG','MF500','MF850','MR2AVG','MR500', 'MR850','PWAVG','RH2AVG','RH500','RH850','SLPAVG','T2AVG','T2MAX','T2MIN','T500','T850','U10AVG', 'U200','U500','U850','V10AVG','V200','V500','V850','WSPD10','WSPD200','WSPD500','WSPD850','Z500']
    rgsWTvarsAll   =VarsFullNameAll
    iSelVariables=[VarsFullNameAll.index(VarsFullName[ii]) for ii in range(len(VarsFullName))]
    rgsWTvars=np.array(rgsWTvarsAll)[np.array(iSelVariables).astype('int')]
    rgsWTfolders=['/glade/campaign/mmm/c3we/prein/ERA5/'+str(VarsFullName[va])+'/' for va in range(len(VarsFullName))]
    # rgsWTfolders=np.array(rgsWTfoldersAll)[np.array(iSelVariables).astype('int')]

    DomDegreeAdd=np.array([2,   5,  10, 20])[['S','M','L','XXL'].index(WT_Domains)]


    return rgdTime, iMonths, sPlotDir, sDataDir, sSubregionPR, rgsWTvars, VarsFullName,rgsWTfolders, rgrNrOfExtremes, WT_Domains, DomDegreeAdd, Annual_Cycle, SpatialSmoothing, Metrics, Dimensions, FireObs, REANAL, ClusterMeth, ClusterBreakup, RelAnnom, NormalizeData, MinDistDD, RemoveAnnualCycl
