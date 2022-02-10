# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 04:59:24 2020


CHIRPS-GEFS [ historical processing]
@author: laver1
"""

import os, sys, os.path
from osgeo import gdal
from geotiff import GeoTiff

import numpy as np
import pandas as pd
import Processtools as process
from skimage.transform import resize
import Verification as vr

from ftplib import FTP
import netCDF4 as nc
from shutil import copyfile
import pickle as pkl

CHIRPSdir = r'E:\CHIRPS-GEFS\Updated_chirps-gefs'
Outdir = r'E:\CHIRPS-GEFS\Process_updated'

geo_tiff = GeoTiff(r'D:\ADPC\landslides\chirps_GEFS_analysis\mask_LMB_01.tif')
Mask = np.array(geo_tiff.read())
mask = Mask == 1
Mask[Mask == 0] = np.nan
#(93,110,9,23)
#%% ###########################################################################
# #######              PROCESS CHIRPS GEFS         ########################### 
# ############################################################################

DateStart='2016-05-01'
DateEnd = '2019-12-31'
Dates_ =  pd.date_range(DateStart,DateEnd, freq='d')

Dates = Dates_[(Dates_.month >= 6) & (Dates_.month <= 10)] 

for Date in Dates:
    
    # Date = Dates[0]    
    folder_path = '{0:02d}/{1:02d}/{2:02d}'.format(Date.year,Date.month,Date.day)        
    # Files = ftp.nlst('cgefs_precip_a0p05/daily_new/archived/'+ folder_path)     
    Files = os.listdir(os.path.join(CHIRPSdir,'{0:02d}\\{1:02d}\\{2:02d}'.format(Date.year,Date.month,Date.day)))

    Day = 1
    for ftpPath_ in Files:
        # ftpPath_ = Files[0]        
        # [Fold,filename] = ftpPath_.split(folder_path+'/')  
        filename = ftpPath_
        MK_files = os.path.join(Outdir,'D{}'.format(Day))      
        chirpsgefs_file = os.path.join(MK_files,filename)
        
        src = os.path.join(CHIRPSdir,'{0:02d}\\{1:02d}\\{2:02d}'.format(Date.year,Date.month,Date.day),filename)
        copyfile(src, chirpsgefs_file)
        Day = Day +1 
        print(Date)              

#%% ###########################################################################
# #######              TEMPORAL ANALYSIS        ########################### 
# ############################################################################

CHIRPSdir = Outdir
IMERG = r'E:\GPM\bico_imerg'

Perfom_chirps = {}

for forecast in range(1,6):
    # forecast = 5
    Files_  = [f for f in os.listdir(os.path.join(CHIRPSdir,'D'+ str(forecast))) if f.endswith( '.nc' )]    
    Mons_monsoon =['.01','.02','.03','.04','.05','.11','.12']
    Files = [x for x in Files_ if
              all(y not in x for y in Mons_monsoon)] # select just in monsson 
    
    CHIRPS = []
    
    for f in Files:
        [Year,MD,ex] = f.replace('mb_cgefs_precip_0p05_','').split('.') 
#        f= Files[0]        
        # CHIRPS-GEFS
        try:                          
                        
            t_imerg = nc.Dataset(os.path.join(IMERG,Year,'rain_'+ Year + MD + '.nc')) # read file IMERG
            imerg = t_imerg.variables['precip'] 
            ILat_min, ILat_max  = np.round( [t_imerg.variables['lat'][:].min(),t_imerg.variables['lat'][:].max()],2)
            ILont_min, ILon_max  = np.round([t_imerg.variables['lon'][:].min(),t_imerg.variables['lon'][:].max()],2)            
            imerg = np.flipud(imerg) *1.14 # transform to day           
            imerg = imerg[0,:,:]
            t_imerg.close()            
            
            t_chirps = nc.Dataset(os.path.join(CHIRPSdir,'D' + str(forecast),f)) # read file
            Lat = t_chirps.variables['lat'][:]
            Lon = t_chirps.variables['lon'][:]
            Lon_index = np.asanyarray(np.where( (Lon >= ILont_min ) & (Lon <= ILon_max) ) )
            Lat_index = np.asanyarray(np.where( (Lat >= ILat_min ) & (Lat <= ILat_max) ) )
            chirps = t_chirps.variables['precip'][Lat_index.min():Lat_index.max()+1, Lon_index.min():Lon_index.max()+1]
            # chirps = t_chirps.variables['precip']
            chirps = np.flipud(chirps)
            t_chirps.close()                 
                
            chirps_scaled = resize(chirps, imerg.shape)
            # imerg_scaled = imerg     
            
            # import matplotlib.pyplot as plt
            # plt.figure(num=None, figsize=(4, 4), dpi=100, facecolor='w', edgecolor='k')
            # plt.subplot(121)
            # plt.imshow(chirps_scaled,vmin=0,vmax=25)
            # plt.subplot(122)
            # plt.imshow(imerg,vmin=0,vmax=25)                             
            
            # error chirps
            chirps_rmse  = np.round(vr.Error.RMSE(imerg.ravel(),chirps_scaled.ravel()),2)
            chirps_bias  =  np.round(vr.Error.BIAS(imerg.ravel(),chirps_scaled.ravel()),2)
            chirps_cc = np.round(np.corrcoef(imerg.ravel(),chirps_scaled.ravel())[1,0],2)
            
            FBI_chirps,POD_chirps,FAR_chirps,POFD_chirps,CSI_chirps  = np.round(vr.Error.cat_linear(imerg.ravel(),chirps_scaled.ravel(),1),2)             
            chirps_error = np.array([chirps_rmse, chirps_bias,chirps_cc,FBI_chirps, POD_chirps, FAR_chirps,CSI_chirps,int(Year)])      
            CHIRPS.append(chirps_error)    
            print(f)
            
        except:
            print('no file' + f)
            
    # errors 
    CHIRPS_list = pd.DataFrame(np.asarray(CHIRPS),columns = ['rmse', 'bias','cc','FBI', 'POD', 'FAR','CSI','Year'])  
    CHIRPS_list.to_csv( os.path.join(r'E:\CHIRPS-GEFS\Performance_analysis','CHIRPS_{}MRC.csv'.format(forecast)))       
    
    Perfom_chirps[forecast] = CHIRPS_list
        
output = open(os.path.join(r'E:\CHIRPS-GEFS\Performance_analysis','Performance_summary_.pkl'), 'wb')
pkl.dump(Perfom_chirps, output)
output.close()            
#%%
#%%  PLOT  ERRORS    
    
forecast = []
Type = []
YY= []
RMSEc,BIASc,FBIc,PODc,FARc,CSIc, CCc = [[],[],[],[],[],[],[]]

for  day in range(2,6): 
#    day = 1
    for val in range(len(Perfom_chirps[day])):
        
        RMSEc.append(Perfom_chirps[day].rmse[val])        
        BIASc.append(Perfom_chirps[day].bias[val])
        CCc.append(Perfom_chirps[day].cc[val])
        FBIc.append(Perfom_chirps[day].FBI[val])
        PODc.append(Perfom_chirps[day].POD[val])
        FARc.append(Perfom_chirps[day].FAR[val])
        CSIc.append(Perfom_chirps[day].CSI[val])
        YY.append(Perfom_chirps[day].Year[val])    
        Type.append('CHIRPS-GEFS')        
        forecast.append(day-1)            
    
   
PERFORM = pd.DataFrame ({'Dataset': np.asarray(Type),'forecast': np.asarray(forecast),
           'RMSE': np.asarray(RMSEc),'BIAS': np.asarray(BIASc),
           'FBI': np.asarray(FBIc),'POD': np.asarray(PODc),
           '1-FAR': 1-np.asarray(FARc),'CSI': np.asarray(CSIc),'CC': np.asarray(CCc),'Year': np.asarray( YY) } )
    
import seaborn as sns
sns.set()
from matplotlib import style
style.use('ggplot')

palette = sns.color_palette("mako_r", 6)

ax = sns.lineplot(x="forecast", y="RMSE", hue="Year", style="Year", markers=True, data=PERFORM,palette=palette)  
ax.set(ylim=(10, None))   

ax = sns.lineplot(x="forecast", y="BIAS", hue="Year", style="Year", markers=True, data=PERFORM,palette=palette)  
 
ax = sns.lineplot(x="forecast", y="FBI", hue="Year", style="Year", markers=True, data=PERFORM,palette=palette) 

ax = sns.lineplot(x="forecast", y="POD", hue="Year", style="Year", markers=True, data=PERFORM,palette=palette)  
ax.set(ylim=(0.5, 1))   

ax = sns.lineplot(x="forecast", y="1-FAR", hue="Year", style="Year", markers=True, data=PERFORM,palette=palette)    
ax.set(ylim=(0.5, 1)) 

ax = sns.lineplot(x="forecast", y="CSI", hue="Year", style="Year", markers=True, data=PERFORM,palette=palette)   
ax.set(ylim=(0.5, 1)) 

ax = sns.lineplot(x="forecast", y="CC", hue="Year", style="Year", markers=True, data=PERFORM,palette=palette)  
ax.set(ylim=(0.2, 0.5)) 




###############################################################################
#%%           2.    INTENSITY ANALYSIS
###############################################################################

perform_output = r'E:\CHIRPS-GEFS\Performance_analysis\intensity'

Intensity = [0.5,1.0,2.5,5,7.5,10,12,15,20,25,30,40]

Perfom_chirps_FBI, Perfom_chirps_POD, Perfom_chirps_FAR, Perfom_chirps_CSI = [{},{},{},{}]   


for forecast in range(1,6):
    # forecast = 1
    Files_  = [f for f in os.listdir(os.path.join(CHIRPSdir,'D'+ str(forecast))) if f.endswith( '.nc' )]    
    Mons_monsoon =['.01','.02','.03','.04','.05','.11','.12']
    Files = [x for x in Files_ if
              all(y not in x for y in Mons_monsoon)] # select just in monsson 
    
    FBI_CH,POD_CH, FAR_CH,CSI_CH,YEAR_CH = [[],[],[],[],[]]          

    
    for f in Files:
        [Year,MD,ex] = f.replace('mb_cgefs_precip_0p05_','').split('.') 
 #        f= Files[0]        
        # CHIRPS-GEFS
        try:                
            t_imerg = nc.Dataset(os.path.join(IMERG,Year,'rain_'+ Year + MD + '.nc')) # read file IMERG
            imerg = t_imerg.variables['precip'] 
            ILat_min, ILat_max  = np.round( [t_imerg.variables['lat'][:].min(),t_imerg.variables['lat'][:].max()],2)
            ILont_min, ILon_max  = np.round([t_imerg.variables['lon'][:].min(),t_imerg.variables['lon'][:].max()],2)            
            imerg = np.flipud(imerg) *1.14 # transform to day           
            imerg = imerg[0,:,:]
            t_imerg.close()            
            
            t_chirps = nc.Dataset(os.path.join(CHIRPSdir,'D' + str(forecast),f)) # read file
            Lat = t_chirps.variables['lat'][:]
            Lon = t_chirps.variables['lon'][:]
            Lon_index = np.asanyarray(np.where( (Lon >= ILont_min ) & (Lon <= ILon_max) ) )
            Lat_index = np.asanyarray(np.where( (Lat >= ILat_min ) & (Lat <= ILat_max) ) )
            chirps = t_chirps.variables['precip'][Lat_index.min():Lat_index.max()+1, Lon_index.min():Lon_index.max()+1]
            # chirps = t_chirps.variables['precip']
            chirps = np.flipud(chirps)
            t_chirps.close()                 
            
                
            chirps_scaled = resize(chirps, imerg.shape)
            imerg_scaled_V = imerg[mask]
            chirps_scaled_V = chirps_scaled[mask]        
            
            # import matplotlib.pyplot as plt
            # plt.figure(num=None, figsize=(4, 4), dpi=100, facecolor='w', edgecolor='k')
            # plt.subplot(131)
            # plt.imshow(chirps,vmin=0,vmax=25,cmap='jet')
            # plt.subplot(132)
            # plt.imshow(imerg,vmin=0,vmax=25,cmap='jet')     
            # plt.subplot(133)
            # plt.imshow(mask)     
            
 
            chirps_error = []
            
            for inten in Intensity:
            
                # error chirps
                FBI_chirps,POD_chirps,FAR_chirps,POFD_chirps,CSI_chirps  = np.round(vr.Error.cat_linear(imerg_scaled_V.ravel(),chirps_scaled_V.ravel(),inten),2)        
                chirps_error_ = np.array([FBI_chirps, POD_chirps, 1-FAR_chirps,CSI_chirps,int(Year)])                     
                
                chirps_error.append(chirps_error_)                 
                
            chirps_error = pd.DataFrame(np.asarray(chirps_error),columns = ['FBI', 'POD', 'FAR','CSI','Year'])               
            
            FBI_CH.append(chirps_error.FBI) 
            POD_CH.append(chirps_error.POD)   
            FAR_CH.append(chirps_error.FAR)   
            CSI_CH.append(chirps_error.CSI)            
            YEAR_CH.append(chirps_error.Year)  
           
            print(f)            
        except:
            print('no file' + f)   
            
    CHIRPS_YEAR =  pd.DataFrame(pd.concat(YEAR_CH, axis=1).values.T[:,0],columns=['year'])    
    
    CHIRPS_FBI = pd.DataFrame(pd.concat(FBI_CH, axis=1).values.T,columns =Intensity)
    CHIRPS_FBI = pd.concat([CHIRPS_FBI,CHIRPS_YEAR], axis=1)
    
    CHIRPS_POD = pd.DataFrame(pd.concat(POD_CH, axis=1).values.T,columns =Intensity) 
    CHIRPS_POD = pd.concat([CHIRPS_POD,CHIRPS_YEAR], axis=1)
    
    CHIRPS_FAR = pd.DataFrame(pd.concat(FAR_CH, axis=1).values.T,columns =Intensity)  
    CHIRPS_FAR = pd.concat([CHIRPS_FAR,CHIRPS_YEAR], axis=1)
    
    CHIRPS_CSI = pd.DataFrame(pd.concat(CSI_CH, axis=1).values.T,columns =Intensity)   
    CHIRPS_CSI = pd.concat([CHIRPS_CSI,CHIRPS_YEAR], axis=1)           
            
    # errors   
    CHIRPS_FBI.to_csv( os.path.join(perform_output,'CHIRPS_FBI_{0}.csv'.format(forecast))) 
    CHIRPS_POD.to_csv( os.path.join(perform_output,'CHIRPS_POD_{0}.csv'.format(forecast))) 
    CHIRPS_FAR.to_csv( os.path.join(perform_output,'CHIRPS_FAR_{0}.csv'.format(forecast))) 
    CHIRPS_CSI.to_csv( os.path.join(perform_output,'CHIRPS_CSI_{0}.csv'.format(forecast)))        

    
    Perfom_chirps_FBI[forecast] = CHIRPS_FBI
    Perfom_chirps_POD[forecast] = CHIRPS_POD
    Perfom_chirps_FAR[forecast] = CHIRPS_FAR
    Perfom_chirps_CSI[forecast] = CHIRPS_CSI    
    

FBIc,PODc,FARc,CSIc = [[],[],[],[]]
for  day in range(1,6):      
    FBIc.append(Perfom_chirps_FBI[day].median())  
    PODc.append(Perfom_chirps_POD[day].median())
    FARc.append(Perfom_chirps_FAR[day].median())
    CSIc.append(Perfom_chirps_CSI[day].median())       
FBIc = pd.concat(FBIc,axis=1) 
PODc = pd.concat(PODc,axis=1)
FARc = pd.concat(FARc,axis=1)
CSIc = pd.concat(CSIc,axis=1)  

#%%      PLOT CHIRPS  

import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
    
from mpl_toolkits.axes_grid1 import ImageGrid


FARc = FARc.drop(index='year')
PODc = PODc.drop(index='year')
CSIc = CSIc.drop(index='year') 
        
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

Xtick = [0,1,2,3]
Ytick = np.linspace(0,len(Intensity),4).astype(int)
Ylabel = np.linspace(min(Intensity),max(Intensity),4).astype(int)
Xlabel = [1,2,3,4]

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


FAR_v = FARc.values
FAR_v[FAR_v==1] = 0

FBI_v = FBIc.values.T
FBI_v[FBI_v==1] = 0

POD_v =PODc.values.T
POD_v[POD_v==1] = 0

CSI_v =CSIc.values.T
CSI_v[FAR_v==1] = 0


fig= plt.figure(num=None, figsize=plt.figaspect(1.0), dpi=150, facecolor='w', edgecolor='k')
## FBI
ax_0 = fig.add_subplot(4, 1, 1)
po_0 = ax_0.imshow(FBIc.values.T, vmin=0, vmax=1.5,cmap = "jet",origin='lower')

##   POD
ax_2 = fig.add_subplot(4, 1, 2)
po_2 = ax_2.imshow(PODc.values.T, vmin=0, vmax=1,cmap = "jet",origin='lower')

##   FAR
ax_4 = fig.add_subplot(4, 1, 3)
po_4 = ax_4.imshow(FARc.values.T, vmin=0, vmax=1,cmap = "jet",origin='lower')

##   CSI
ax_6 = fig.add_subplot(4, 1, 4)
po_6 = ax_6.imshow(CSIc.values.T, vmin=0, vmax=1,cmap = "jet",origin='lower')


ax_0.set_title('CHIRPS-GEFS')
ax_0.set_ylabel('FBI')
ax_2.set_ylabel('POD')
ax_4.set_ylabel('FAR')
ax_6.set_ylabel('CSI')
ax_6.set_xlabel('intensity ($mm$)')


Ytick = [0,2,4,5]
Xtick = [0,5,10,15]
Xlabel = [0.5,10,40,80]

ax_0.set_yticks(Ytick )
ax_0.set_xticks(Xtick)
ax_2.set_yticks(Ytick)
ax_2.set_xticks(Xtick)
ax_4.set_yticks(Ytick)
ax_4.set_xticks(Xtick)
ax_6.set_yticks(Ytick)
ax_6.set_xticks(Xtick)
ax_6.set_xticklabels(Xlabel)
ax_0.set_xticklabels([ ])
ax_2.set_xticklabels([ ])
ax_4.set_xticklabels([ ])


###############################################################################
#%%           3.    SPATIAL ANALSYS
###############################################################################

perform_output = r'E:\CHIRPS-GEFS\Performance_analysis\Spatial'
inten = 1

Perfom_chirps_FBI, Perfom_chirps_POD, Perfom_chirps_FAR, Perfom_chirps_CSI = [{},{},{},{}]   

for forecast in range(1,6):
    # forecast = 1
    Files_  = [f for f in os.listdir(os.path.join(CHIRPSdir,'D'+ str(forecast))) if f.endswith( '.nc' )]    
    Mons_monsoon =['.01','.02','.03','.04','.05','.11','.12']
    Files = [x for x in Files_ if
              all(y not in x for y in Mons_monsoon)] # select just in monsson 
    
    [matrix_chirps, matrix_gfs, matrix_imerg]= [[],[],[]]     

    
    for f in Files:
        [Year,MD,ex] = f.replace('mb_cgefs_precip_0p05_','').split('.') 
 #        f= Files[0]        
        # CHIRPS-GEFS
        try:  
            t_imerg = nc.Dataset(os.path.join(IMERG,Year,'rain_'+ Year + MD + '.nc')) # read file IMERG
            imerg = t_imerg.variables['precip'] 
            ILat_min, ILat_max  = np.round( [t_imerg.variables['lat'][:].min(),t_imerg.variables['lat'][:].max()],2)
            ILont_min, ILon_max  = np.round([t_imerg.variables['lon'][:].min(),t_imerg.variables['lon'][:].max()],2)            
            imerg = np.flipud(imerg) *1.14 # transform to day           
            imerg = imerg[0,:,:]
            t_imerg.close()            
            
            t_chirps = nc.Dataset(os.path.join(CHIRPSdir,'D' + str(forecast),f)) # read file
            Lat = t_chirps.variables['lat'][:]
            Lon = t_chirps.variables['lon'][:]
            Lon_index = np.asanyarray(np.where( (Lon >= ILont_min ) & (Lon <= ILon_max) ) )
            Lat_index = np.asanyarray(np.where( (Lat >= ILat_min ) & (Lat <= ILat_max) ) )
            chirps = t_chirps.variables['precip'][Lat_index.min():Lat_index.max()+1, Lon_index.min():Lon_index.max()+1]
            # chirps = t_chirps.variables['precip']
            chirps = np.flipud(chirps)
            t_chirps.close()     
    
            chirps_scaled = resize(chirps, imerg.shape)
            
            matrix_chirps.append(chirps_scaled)
            matrix_imerg.append(imerg)
            print(f)            
        except:
            print('no file' + f)   
            
    matrix_chirps = np.array(matrix_chirps)   
    matrix_imerg  = np.array(matrix_imerg) 
    
    # categorical spatial error
    [T,X,Y] = matrix_chirps.shape
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

#%% PLOT
import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
    
from mpl_toolkits.axes_grid1 import ImageGrid
import shapefile

extent = (ILont_min, ILon_max, ILat_min,ILat_max) # extent of the plot
#  Area 
sf = shapefile.Reader(r'D:\ADPC\Data\Shape_Mekong\Mekong_basin_area.shp')
shapesB = sf.shape(0)
intervalsB = list(shapesB.parts) + [len(shapesB.points)]
pointsB = np.array(shapesB.points)

sf = shapefile.Reader(r'D:\ADPC\Data\Shape_Mekong\Mekong_land_areaSimpl.shp')
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
    
    grid = ImageGrid(fig, 111,  nrows_ncols=(2,2), axes_pad=0.35,  share_all=True,
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
    
    grid[3].cax.colorbar(po)
    grid[3].cax.toggle_label(True)
    grid[0].set_ylabel('lat')
    grid[2].set_ylabel('lat')
    grid[2].set_xlabel('lon')
    grid[3].set_xlabel('lon')

spa_plot(Perfom_chirps_POD)
spa_plot(Perfom_chirps_FAR)
spa_plot(Perfom_chirps_CSI)    


#%%
 