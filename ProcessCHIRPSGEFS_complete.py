# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 04:59:24 2020


CHIRPS-GEFS [ historical processing]
@author: laver1
"""

import os, sys, os.path
from osgeo import gdal
import numpy as np
import pandas as pd
import Processtools as process
from skimage.transform import resize
import Verification as vr

RasterFormat = 'GTiff'
PixelRes = 0.05 #degrees

MinLon = 93
MaxLon = 110
MinLat = 9
MaxLat = 23 # 35

EPSG='EPSG:4326'
BoundariesP = [MinLon,MinLat,MaxLon,MaxLat]

Mask = gdal.Open(r'D:\PhD\ADPC\Data\Raster_LMB\MASK_LMB_0.25.tif').ReadAsArray()
mask = Mask == 1
Mask[Mask == 0] = np.nan
#(93,110,9,23)
#%% ###########################################################################
# #######              PROCESS CHIRPS GEFS         ########################### 
# ############################################################################

CHIRPSdir = r'G:\CHIRPS-GEFS\last'
Outdir = r'G:\CHIRPS-GEFS\Process_Last'

for year in ['2017','2018','2019']:
    Dates =   [f[5:14] for f in os.listdir(os.path.join(CHIRPSdir,year))]
    Creation = [f[-13:-4] for f in os.listdir(os.path.join(CHIRPSdir,year))]
    
    Date_creation = np.unique(np.asarray(Creation).astype(float))
    
    for no in range(len(Date_creation)):
    #    no = 0    
        File  = [f for f in os.listdir(os.path.join(CHIRPSdir,year)) if f.endswith( '{:05}.tif'.format(Date_creation[no]) )]
        
        for day in range(len(File)):
    #        day = 0
            Raster = gdal.Open(os.path.join(CHIRPSdir,year,File[day]))
            Raster_save = gdal.Warp(os.path.join(Outdir,str(day),File[day]),Raster, xRes=PixelRes, yRes=PixelRes,outputBounds=BoundariesP)
            print(File[day])
#Raster3= gdal.Open(r'E:\CHIRPS-GEFS\\2019\data.2019.0614.created-from.2019.0611.tif')

#%% ###########################################################################
# #######              TEMPORAL ANALYSIS        ########################### 
# ############################################################################
            
CHIRPSdir = r'G:\CHIRPS-GEFS\Process_Last'
GFSdir = r'G:\GFS'        
gfs_for = ['a','b','c','d','e','f','g']
IMERG = r'G:\E_IMERG_DT\BIASCor_FINAL'
perform_output = r'G:\PerfomarceOutput\temporal'
Perfom_chirps = {}
Perform_gfs = {}

for forecast in range(10):    
    print(str(forecast))
        
    File  = [f for f in os.listdir(os.path.join(CHIRPSdir,str(forecast))) if f.endswith( '.tif' )]
    
    CHIRPS = []
    GFS = []
    for f in File:
#        Date = f[10:14] 
#        f= File[18]        
        year = f[7:9] 
        # IMERG
        imerg = gdal.Open(os.path.join(IMERG,'rain_20'+ year + f[10:14] +  '000000.asc')).ReadAsArray()*1.14 # transform to day
        imerg_scaled = resize(imerg, [56,68]) #gfs.shape    
        imerg_scaled_V = imerg_scaled[mask]
        
        try:
            # CHIRPS-GEFS
            chirps = gdal.Open(os.path.join(CHIRPSdir,str(forecast),f)).ReadAsArray()     
            chirps_scaled = resize(chirps, [56,68]) #gfs.shape  
            chirps_scaled_V = chirps_scaled[mask]
              
            # error chirps
            chirps_rmse  = np.round(vr.Error.RMSE(imerg_scaled_V.ravel(),chirps_scaled_V.ravel()),2)
            chirps_bias  =  np.round(vr.Error.BIAS(imerg_scaled_V.ravel(),chirps_scaled_V.ravel()),2)
            chirps_cc    =  np.round(np.corrcoef(imerg_scaled_V.ravel(),chirps_scaled_V.ravel()),2)[1,0]
            FBI_chirps,POD_chirps,FAR_chirps,POFD_chirps,CSI_chirps  = np.round(vr.Error.cat_linear(imerg_scaled_V.ravel(),chirps_scaled_V.ravel(),1),2) 
            
            chirps_error = np.array([chirps_rmse, chirps_bias,FBI_chirps, POD_chirps, FAR_chirps,CSI_chirps, chirps_cc]) 
            CHIRPS.append(chirps_error)              
        except:
            chirps_error = np.array([np.nan, np.nan, np.nan, np.nan, np.nan,np.nan, np.nan])          
            CHIRPS.append(chirps_error)
            print('no file chirps' + f)
            
        try:            
            # GFS
            gfs = gdal.Open(os.path.join(GFSdir,year,year+'_'+ gfs_for[forecast],'_' + year + f[10:14] + gfs_for[forecast] +  '.tif')).ReadAsArray()
            gfs_V = gfs[mask]
            
            # error gfs
            gfs_rmse  = np.round(vr.Error.RMSE(imerg_scaled_V.ravel(),gfs_V.ravel()),2)
            gfs_bias  =  np.round(vr.Error.BIAS(imerg_scaled_V.ravel(),gfs_V.ravel()),2)
            gfs_cc    =  np.round(np.corrcoef(imerg_scaled_V.ravel(),gfs_V.ravel()),2)[1,0] 
            FBI_gfs,POD_gfs,FAR_gfs,POFD_gfs,CSI_gfs  = np.round(vr.Error.cat_linear(imerg_scaled_V.ravel(),gfs_V.ravel(),1),2) 
    
            gfs_error = np.array([gfs_rmse, gfs_bias, FBI_gfs, POD_gfs, FAR_gfs,CSI_gfs, gfs_cc])          
            GFS.append(gfs_error)            
            print(f)            
        except:
            gfs_error = np.array([np.nan, np.nan, np.nan, np.nan, np.nan,np.nan, np.nan])          
            GFS.append(gfs_error)
            print('no file_gfs' + f)
            
    # errors 
    CHIRPS_list = pd.DataFrame(np.asarray(CHIRPS),columns = ['rmse', 'bias','FBI', 'POD', 'FAR','CSI','cc'])  
    CHIRPS_list.to_csv( os.path.join(perform_output,'CHIRPS_{}MRC.csv'.format(forecast)))        
    
    GFS_list = pd.DataFrame(np.asarray(GFS),columns = ['rmse', 'bias','FBI', 'POD', 'FAR','CSI','cc'])  
    GFS_list.to_csv( os.path.join(perform_output,'GFS_{}MRC.csv'.format(forecast)))  
    
    Perfom_chirps[forecast] = CHIRPS_list
    Perform_gfs[forecast] =  GFS_list          
            
#%%  PLOT  ERRORS    
    
forecast = []
Type = []
RMSEc,BIASc,FBIc,PODc,FARc,CSIc, CCc = [[],[],[],[],[],[],[]]

for  day in range(1,10): 
#    day = 1
    for val in range(len(Perfom_chirps[day])):
        
        RMSEc.append(Perfom_chirps[day].rmse[val])        
        BIASc.append(Perfom_chirps[day].bias[val])
        CCc.append(Perfom_chirps[day].cc[val])
        FBIc.append(Perfom_chirps[day].FBI[val])
        PODc.append(Perfom_chirps[day].POD[val])
        FARc.append(Perfom_chirps[day].FAR[val])
        CSIc.append(Perfom_chirps[day].CSI[val])
        Type.append('CHIRPS-GEFS')
        forecast.append(day)
        
        RMSEc.append(Perform_gfs[day].rmse[val])
        BIASc.append(Perform_gfs[day].bias[val])
        CCc.append(Perform_gfs[day].cc[val])
        FBIc.append(Perform_gfs[day].FBI[val])
        PODc.append(Perform_gfs[day].POD[val])
        FARc.append(Perform_gfs[day].FAR[val])
        CSIc.append(Perform_gfs[day].CSI[val])
        Type.append('GFS')        
        forecast.append(day)        
    
PERFORM = pd.DataFrame ({'Dataset': np.asarray(Type),'forecast': np.asarray(forecast),
           'RMSE': np.asarray(RMSEc),'BIAS': np.asarray(BIASc),
           'FBI': np.asarray(FBIc),'POD': np.asarray(PODc),
           '1-FAR': 1-np.asarray(FARc),'CSI': np.asarray(CSIc),'CC': np.asarray(CCc) } )
    
import seaborn as sns
sns.set()
from matplotlib import style
style.use('ggplot')

ax = sns.lineplot(x="forecast", y="RMSE", hue="Dataset", style="Dataset", markers=True, data=PERFORM)    

ax = sns.lineplot(x="forecast", y="BIAS", hue="Dataset", style="Dataset", markers=True, data=PERFORM) 
 
ax = sns.lineplot(x="forecast", y="FBI", hue="Dataset", style="Dataset", markers=True, data=PERFORM)  

ax = sns.lineplot(x="forecast", y="POD", hue="Dataset", style="Dataset", markers=True, data=PERFORM)  

ax = sns.lineplot(x="forecast", y="1-FAR", hue="Dataset", style="Dataset", markers=True, data=PERFORM)  

ax = sns.lineplot(x="forecast", y="CSI", hue="Dataset", style="Dataset", markers=True, data=PERFORM) 

ax = sns.lineplot(x="forecast", y="CC", hue="Dataset", style="Dataset", markers=True, data=PERFORM) 

###############################################################################
#%%           2.    INTENSITY ANALYSIS
###############################################################################

CHIRPSdir = r'G:\CHIRPS-GEFS\Process_Last'
GFSdir = r'G:\GFS'        
gfs_for = ['a','b','c','d','e','f','g']
IMERG = r'G:\E_IMERG_DT\BIASCor_FINAL'
perform_output = r'G:\PerfomarceOutput\Intensity'

Intensity = [0.5,1.0,2.5,5,7.5,10,12,15,20,25,30,40]

Perfom_chirps_FBI, Perfom_chirps_POD, Perfom_chirps_FAR, Perfom_chirps_CSI = [{},{},{},{}]   
Perform_gfs_FBI, Perform_gfs_POD, Perform_gfs_FAR, Perform_gfs_CSI = [{},{},{},{}]   

for forecast in range(7):
    
    File  = [f for f in os.listdir(os.path.join(CHIRPSdir,str(forecast))) if f.endswith( '.tif' )]
    
    FBI_CH,POD_CH, FAR_CH,CSI_CH = [[],[],[],[]]            
    FBI_GS,POD_GS,FAR_GS,CSI_GS = [[],[],[],[]]   
    for f in File:
#        Date = f[10:14] 
#        f= File[0] 
        year = f[7:9]         
        try:
            # CHIRPS-GEFS
            chirps = gdal.Open(os.path.join(CHIRPSdir,str(forecast),f)).ReadAsArray()            
            # GFS
            gfs = gdal.Open(os.path.join(GFSdir,year,year+'_'+ gfs_for[forecast],'_' + year + f[10:14] + gfs_for[forecast] +  '.tif')).ReadAsArray()
            
            imerg = gdal.Open(os.path.join(IMERG,'rain_20'+ year + f[10:14] +  '000000.asc')).ReadAsArray()*1.14 # transform to day
    
            chirps_scaled = resize(chirps, gfs.shape)
            imerg_scaled = resize(imerg, gfs.shape)
            
            imerg_scaled_V = imerg_scaled[mask]
            chirps_scaled_V = chirps_scaled[mask]
            gfs_V = gfs[mask]
 
            chirps_error = []
            gfs_error = []
            for inten in Intensity:
            
                # error chirps
                FBI_chirps,POD_chirps,FAR_chirps,POFD_chirps,CSI_chirps  = np.round(vr.Error.cat_linear(imerg_scaled_V.ravel(),chirps_scaled_V.ravel(),inten),2)        
                chirps_error_ = np.array([FBI_chirps, POD_chirps, 1-FAR_chirps,CSI_chirps])           
    
                # error gfs
                FBI_gfs,POD_gfs,FAR_gfs,POFD_gfs,CSI_gfs  = np.round(vr.Error.cat_linear(imerg_scaled_V.ravel(),gfs_V.ravel(),inten),2)     
                gfs_error_ = np.array([FBI_gfs, POD_gfs, 1-FAR_gfs,CSI_gfs]) 
                
                chirps_error.append(chirps_error_) 
                gfs_error.append(gfs_error_) 
                
            chirps_error = pd.DataFrame(np.asarray(chirps_error),columns = ['FBI', 'POD', 'FAR','CSI']) 
            gfs_error = pd.DataFrame(np.asarray(gfs_error),columns = ['FBI', 'POD', 'FAR','CSI'])  
            
            FBI_CH.append(chirps_error.FBI) 
            POD_CH.append(chirps_error.POD)   
            FAR_CH.append(chirps_error.FAR)   
            CSI_CH.append(chirps_error.CSI)  
            
            FBI_GS.append(gfs_error.FBI) 
            POD_GS.append(gfs_error.POD)   
            FAR_GS.append(gfs_error.FAR)  
            CSI_GS.append(gfs_error.CSI)             
            print(f)            
        except:
            print('no file' + f)   
            
    CHIRPS_FBI = pd.DataFrame(pd.concat(FBI_CH, axis=1).values.T,columns =Intensity)
    CHIRPS_POD = pd.DataFrame(pd.concat(POD_CH, axis=1).values.T,columns =Intensity) 
    CHIRPS_FAR = pd.DataFrame(pd.concat(FAR_CH, axis=1).values.T,columns =Intensity)  
    CHIRPS_CSI = pd.DataFrame(pd.concat(CSI_CH, axis=1).values.T,columns =Intensity)
    
    GFS_FBI = pd.DataFrame(pd.concat(FBI_GS, axis=1).values.T,columns =Intensity)
    GFS_POD = pd.DataFrame(pd.concat(POD_GS, axis=1).values.T,columns =Intensity) 
    GFS_FAR = pd.DataFrame(pd.concat(FAR_GS, axis=1).values.T,columns =Intensity)
    GFS_CSI = pd.DataFrame(pd.concat(CSI_GS, axis=1).values.T,columns =Intensity)               
            
    # errors   
    CHIRPS_FBI.to_csv( os.path.join(perform_output,'CHIRPS_FBI_{0}.csv'.format(forecast))) 
    CHIRPS_POD.to_csv( os.path.join(perform_output,'CHIRPS_POD_{0}.csv'.format(forecast))) 
    CHIRPS_FAR.to_csv( os.path.join(perform_output,'CHIRPS_FAR_{0}.csv'.format(forecast))) 
    CHIRPS_CSI.to_csv( os.path.join(perform_output,'CHIRPS_CSI_{0}.csv'.format(forecast)))    
    
    GFS_FBI.to_csv( os.path.join(perform_output,'GFS_FBI_{0}.csv'.format(forecast))) 
    GFS_POD.to_csv( os.path.join(perform_output,'GFS_POD_{0}.csv'.format(forecast))) 
    GFS_FAR.to_csv( os.path.join(perform_output,'GFS_FAR_{0}.csv'.format(forecast))) 
    GFS_CSI.to_csv( os.path.join(perform_output,'GFS_CSI_{0}.csv'.format(forecast))) 

    
    Perfom_chirps_FBI[forecast] = CHIRPS_FBI
    Perfom_chirps_POD[forecast] = CHIRPS_POD
    Perfom_chirps_FAR[forecast] = CHIRPS_FAR
    Perfom_chirps_CSI[forecast] = CHIRPS_CSI
    
    Perform_gfs_FBI[forecast] =  GFS_FBI         
    Perform_gfs_POD[forecast] =  GFS_POD     
    Perform_gfs_FAR[forecast] =  GFS_FAR 
    Perform_gfs_CSI[forecast] =  GFS_CSI      
    

FBIc,PODc,FARc,CSIc = [[],[],[],[]]
FBIg,PODg,FARg,CSIg = [[],[],[],[]]

for  day in range(1,7):  
    
    FBIc.append(Perfom_chirps_FBI[day].median())  
    PODc.append(Perfom_chirps_POD[day].median())
    FARc.append(Perfom_chirps_FAR[day].median())
    CSIc.append(Perfom_chirps_CSI[day].median())    

    FBIg.append(Perform_gfs_FBI[day].median())
    PODg.append(Perform_gfs_POD[day].median())
    FARg.append(Perform_gfs_FAR[day].median())
    CSIg.append(Perform_gfs_CSI[day].median())    
    
FBIc = pd.concat(FBIc,axis=1) 
PODc = pd.concat(PODc,axis=1)
FARc = pd.concat(FARc,axis=1)
CSIc = pd.concat(CSIc,axis=1)  

FBIg = pd.concat(FBIg,axis=1)
PODg = pd.concat(PODg,axis=1)
FARg = pd.concat(FARg,axis=1)
CSIg = pd.concat(CSIg,axis=1) 

#%%      PLOT CHIRPS  

import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
    
from mpl_toolkits.axes_grid1 import ImageGrid


FAR_1 = FARc.values
FAR_1[FAR_1==1] = 0

# Set up figure and image grid
fig = plt.figure(figsize=(9.75, 3))

grid = ImageGrid(fig, 111,  nrows_ncols=(1,3), axes_pad=0.15,  share_all=True,
                 cbar_location="right", cbar_mode="single", cbar_size="7%", cbar_pad=0.15, )

# Add data to image grid
po_2 = grid[0].imshow(PODc.values, vmin=0, vmax=1,origin='lower',cmap = "jet")
po_4 = grid[1].imshow(FAR_1, vmin=0, vmax=1,origin='lower',cmap = "jet")
po_6 = grid[2].imshow(CSIc.values, vmin=0, vmax=1,origin='lower',cmap = "jet")

grid[0].set_ylabel('intensity ($mm$)')
grid[0].set_xlabel('forecast ($days$)')
grid[1].set_xlabel('forecast ($days$)')
grid[2].set_xlabel('forecast ($days$)')
grid[0].set_title('POD')
grid[1].set_title('1-FAR')
grid[2].set_title('CSI')

Xtick = [0,2,4]
Ytick = np.linspace(0,len(Intensity),4).astype(int)
Ylabel = np.linspace(min(Intensity),max(Intensity),4).astype(int)
Xlabel = [1,3,5]

for ax in grid:
    ax.set_xticks(Xtick)
    ax.set_yticks(Ytick)
    ax.set_xticklabels(Xlabel)    
grid[0].set_yticklabels(Ylabel)    

# Colorbar
grid[2].cax.colorbar(po_6)
grid[2].cax.toggle_label(True)

#plt.tight_layout()    # Works, but may still require rect paramater to keep colorbar labels visible
plt.show()

#%%   PLOT GFS

# Set up figure and image grid
fig = plt.figure(figsize=(9.75, 3))

grid = ImageGrid(fig, 111,  nrows_ncols=(1,3), axes_pad=0.15,  share_all=True,
                 cbar_location="right", cbar_mode="single", cbar_size="7%", cbar_pad=0.15, )

# Add data to image grid
po_2 = grid[0].imshow(PODg.values, vmin=0, vmax=1,origin='lower',cmap = "jet")
po_4 = grid[1].imshow(FARg.values, vmin=0, vmax=1,origin='lower',cmap = "jet")
po_6 = grid[2].imshow(CSIg.values, vmin=0, vmax=1,origin='lower',cmap = "jet")
grid[0].set_ylabel('intensity ($mm$)')
grid[0].set_xlabel('forecast ($days$)')
grid[1].set_xlabel('forecast ($days$)')
grid[2].set_xlabel('forecast ($days$)')
grid[0].set_title('POD')
grid[1].set_title('1 - FAR')
grid[2].set_title('CSI')

Xtick = [0,2,4]
Ytick = np.linspace(0,len(Intensity),4).astype(int)
Ylabel = np.linspace(min(Intensity),max(Intensity),4).astype(int)
Xlabel = [1,3,5]

for ax in grid:
    ax.set_xticks(Xtick)
    ax.set_yticks(Ytick)
    ax.set_xticklabels(Xlabel)    
grid[0].set_yticklabels(Ylabel)    

# Colorbar
grid[2].cax.colorbar(po_6)
grid[2].cax.toggle_label(True)

#plt.tight_layout()    # Works, but may still require rect paramater to keep colorbar labels visible
plt.show()

#%%   PLOT COMPLETE
fig= plt.figure(num=None, figsize=plt.figaspect(1.0), dpi=150, facecolor='w', edgecolor='k')
## FBI
ax_0 = fig.add_subplot(4, 2, 1)
po_0 = ax_0.imshow(FBIc.values.T, vmin=0, vmax=1.5)
ax_1 = fig.add_subplot(4, 2, 2)
po_1 = ax_1.imshow(FBIg.values.T, vmin=0, vmax=1.5)

##   POD
ax_2 = fig.add_subplot(4, 2, 3)
po_2 = ax_2.imshow(PODc.values.T, vmin=0, vmax=1)
ax_3 = fig.add_subplot(4, 2, 4)
po_3 = ax_3.imshow(PODg.values.T, vmin=0, vmax=1)

##   FAR
ax_4 = fig.add_subplot(4, 2, 5)
po_4 = ax_4.imshow(FARc.values.T, vmin=0, vmax=1)
ax_5 = fig.add_subplot(4, 2, 6)
po_5 = ax_5.imshow(FARg.values.T, vmin=0, vmax=1)

##   CSI
ax_6 = fig.add_subplot(4, 2, 7)
po_6 = ax_6.imshow(CSIc.values.T, vmin=0, vmax=1)
ax_7 = fig.add_subplot(4, 2, 8)
po_7 = ax_7.imshow(CSIg.values.T, vmin=0, vmax=1)

ax_0.set_title('CHIRPS-GEFS')
ax_1.set_title('GFS')
ax_0.set_ylabel('FBI')
ax_2.set_ylabel('POD')
ax_4.set_ylabel('FAR')
ax_6.set_ylabel('CSI')
ax_6.set_xlabel('intensity ($mm$)')
ax_7.set_xlabel('intensity ($mm$)')

Ytick = [0,2,4,5]
Xtick = [0,5,10,15]
Xlabel = [0.5,10,40,80]

ax_0.set_yticks(Ytick )
ax_0.set_xticks(Xtick)
ax_1.set_yticks(Ytick)
ax_1.set_xticks(Xtick)
ax_2.set_yticks(Ytick)
ax_2.set_xticks(Xtick)
ax_3.set_yticks(Ytick)
ax_3.set_xticks(Xtick)
ax_4.set_yticks(Ytick)
ax_4.set_xticks(Xtick)
ax_5.set_yticks(Ytick)
ax_5.set_xticks(Xtick)
ax_6.set_yticks(Ytick)
ax_6.set_xticks(Xtick)
ax_7.set_yticks(Ytick)
ax_7.set_xticks(Xtick)
ax_6.set_xticklabels(Xlabel)
ax_7.set_xticklabels(Xlabel)
ax_0.set_xticklabels([ ])
ax_1.set_xticklabels([ ])
ax_2.set_xticklabels([ ])
ax_3.set_xticklabels([ ])
ax_4.set_xticklabels([ ])
ax_5.set_xticklabels([ ])
cbar = fig.colorbar(po_1, ax=ax_1, shrink=1,cmap='jet')
cbar = fig.colorbar(po_3, ax=ax_3, shrink=1,cmap='jet')
cbar = fig.colorbar(po_5, ax=ax_5, shrink=1,cmap='jet')
cbar = fig.colorbar(po_7, ax=ax_7, shrink=1,cmap='jet')

###############################################################################
#%%           3.    SPATIAL ANALSYS
###############################################################################

CHIRPSdir = r'G:\CHIRPS-GEFS\Process_Last'
GFSdir = r'G:\GFS'        
gfs_for = ['a','b','c','d','e','f','g']
IMERG = r'G:\E_IMERG_DT\BIASCor_FINAL'
perform_output = r'G:\PerfomarceOutput\spatial'
inten = 1

Perfom_chirps_FBI, Perfom_chirps_POD, Perfom_chirps_FAR, Perfom_chirps_CSI = [{},{},{},{}]   
Perform_gfs_FBI, Perform_gfs_POD, Perform_gfs_FAR, Perform_gfs_CSI = [{},{},{},{}]   

for forecast in range(7):
    
    File  = [f for f in os.listdir(os.path.join(CHIRPSdir,str(forecast))) if f.endswith( '.tif' )]
    [matrix_chirps, matrix_gfs, matrix_imerg]= [[],[],[]]
    
    for f in File:
#        Date = f[10:14] 
#        f= File[0] 
        year = f[7:9]         
        try:
            # CHIRPS-GEFS
            chirps = gdal.Open(os.path.join(CHIRPSdir,str(forecast),f)).ReadAsArray()            
            # GFS
            gfs = gdal.Open(os.path.join(GFSdir,year,year+'_'+ gfs_for[forecast],'_' + year + f[10:14] + gfs_for[forecast] +  '.tif')).ReadAsArray()
            
            imerg = gdal.Open(os.path.join(IMERG,'rain_20'+ year + f[10:14] +  '000000.asc')).ReadAsArray()*1.14 # transform to day
    
            chirps_scaled = resize(chirps, gfs.shape)
            imerg_scaled = resize(imerg, gfs.shape)   
            
            matrix_chirps.append(chirps_scaled)
            matrix_gfs.append(gfs)
            matrix_imerg.append(imerg_scaled)
            print(f)            
        except:
            print('no file' + f)   
            
    matrix_chirps = np.array(matrix_chirps)
    matrix_gfs = np.array(matrix_gfs ) 
    matrix_imerg  = np.array(matrix_imerg) 
    
    # categorical spatial error
    [T,X,Y] = matrix_gfs.shape
    FBI_CH,POD_CH, FAR_CH,CSI_CH = [np.zeros([X,Y]),np.zeros([X,Y]),np.zeros([X,Y]),np.zeros([X,Y])]            
    FBI_GS,POD_GS,FAR_GS,CSI_GS =  [np.zeros([X,Y]),np.zeros([X,Y]),np.zeros([X,Y]),np.zeros([X,Y])]          
        
    for x in range(X):
        for y in range(Y):
            
            # error chirps
            try:                
                FBI_chirps,POD_chirps,FAR_chirps,POFD_chirps,CSI_chirps  = np.round(vr.Error.cat_linear(matrix_imerg[:,x,y],
                                                                                                    matrix_chirps[:,x,y],inten),2)        
                FBI_CH[x,y],POD_CH[x,y], FAR_CH[x,y],CSI_CH[x,y] = [FBI_chirps, POD_chirps,1- FAR_chirps,CSI_chirps]          
            except:
                FBI_CH[x,y],POD_CH[x,y], FAR_CH[x,y],CSI_CH[x,y] = [np.nan, np.nan, np.nan, np.nan]          

            # error gfs
            try:
                FBI_gfs,POD_gfs,FAR_gfs,POFD_gfs,CSI_gfs  = np.round(vr.Error.cat_linear(matrix_imerg[:,x,y],
                                                                                     matrix_gfs[:,x,y],inten),2)     
                FBI_GS[x,y],POD_GS[x,y],FAR_GS[x,y],CSI_GS[x,y] = [FBI_gfs, POD_gfs, 1- FAR_gfs,CSI_gfs]
            except:
                FBI_GS[x,y],POD_GS[x,y],FAR_GS[x,y],CSI_GS[x,y] = [np.nan, np.nan, np.nan, np.nan]
    
    # errors   
    np.savetxt( os.path.join(perform_output,'CHIRPS_FBI_{0}.csv'.format(forecast)), FBI_CH, delimiter=",") 
    np.savetxt( os.path.join(perform_output,'CHIRPS_POD_{0}.csv'.format(forecast)), FBI_CH, delimiter=",") 
    np.savetxt( os.path.join(perform_output,'CHIRPS_FAR_{0}.csv'.format(forecast)), FAR_CH, delimiter=",")
    np.savetxt( os.path.join(perform_output,'CHIRPS_CSI_{0}.csv'.format(forecast)), CSI_CH, delimiter=",")    
    
    np.savetxt( os.path.join(perform_output,'GFS_FBI_{0}.csv'.format(forecast)), CSI_CH, delimiter=",") 
    np.savetxt( os.path.join(perform_output,'GFS_POD_{0}.csv'.format(forecast)), POD_GS, delimiter=",") 
    np.savetxt( os.path.join(perform_output,'GFS_FAR_{0}.csv'.format(forecast)), FAR_GS, delimiter=",") 
    np.savetxt( os.path.join(perform_output,'GFS_CSI_{0}.csv'.format(forecast)), CSI_GS, delimiter=",") 

    Perfom_chirps_FBI[forecast] = FBI_CH
    Perfom_chirps_POD[forecast] = POD_CH
    Perfom_chirps_FAR[forecast] = FAR_CH
    Perfom_chirps_CSI[forecast] = CSI_CH
    
    Perform_gfs_FBI[forecast] =  FBI_GS         
    Perform_gfs_POD[forecast] =  POD_GS     
    Perform_gfs_FAR[forecast] =  FAR_GS 
    Perform_gfs_CSI[forecast] =  CSI_GS   

#%% PLOT
import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
    
from mpl_toolkits.axes_grid1 import ImageGrid
import shapefile

extent = (MinLon, MaxLon, MinLat,MaxLat) # extent of the plot
#  Area 
sf = shapefile.Reader(r'D:\PhD\ADPC\Data\Shape_Mekong\Mekong_basin_area.shp')
shapesB = sf.shape(0)
intervalsB = list(shapesB.parts) + [len(shapesB.points)]
pointsB = np.array(shapesB.points)

sf = shapefile.Reader(r'D:\PhD\ADPC\Data\Shape_Mekong\Mekong_land_areaSimpl.shp')
shapesL = sf.shape(0)
intervalsL = list(shapesL.parts) + [len(shapesL.points)]
pointsL = np.array(shapesL.points)


# interp
from scipy import interpolate
def IntepArray(array):
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    # mask invalid values
    array = np.ma.masked_invalid(array)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]

    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
                          (xx, yy),
                             method='cubic') 
    return GD1

#%%  

def spa_plot(Perfom_chirps_POD):
  
    # Set up figure and image grid
    fig = plt.figure(figsize=(9.75, 3))
    
    grid = ImageGrid(fig, 111,  nrows_ncols=(2,3), axes_pad=0.35,  share_all=True,
                     cbar_location="right", cbar_mode="single", cbar_size="7%", cbar_pad=0.15, )
    
    # Add data to image grid
    for forecast in range(1,len(Perfom_chirps_POD)):        
        po = grid[forecast-1].imshow(IntepArray(Perfom_chirps_POD[forecast])*Mask,extent=extent, 
                 vmin=0, vmax=1,cmap = "jet")#,interpolation='bicubic')
        for (i, j) in zip(intervalsB[:-1], intervalsB[1:]):
             grid[forecast-1].plot(*zip(*pointsB[i:j]),color='grey')
        for (i, j) in zip(intervalsL[:-1], intervalsL[1:]):
             grid[forecast-1].plot(*zip(*pointsL[i:j]),color='black')
    
        grid[forecast-1].set_title('Day_{}'.format(forecast))
    
    grid[5].cax.colorbar(po)
    grid[5].cax.toggle_label(True)
    grid[0].set_ylabel('lat')
    grid[3].set_ylabel('lat')
    grid[3].set_xlabel('lon')
    grid[4].set_xlabel('lon')
    grid[5].set_xlabel('lon')

spa_plot(Perfom_chirps_POD)
spa_plot(Perfom_chirps_FAR)
spa_plot(Perfom_chirps_CSI)    

spa_plot(Perform_gfs_POD)
spa_plot(Perform_gfs_FAR)
spa_plot(Perform_gfs_CSI)

#%%
 