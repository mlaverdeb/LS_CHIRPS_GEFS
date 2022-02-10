# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 04:59:24 2020


CHIRPS-GEFS [ historical processing]
@author: laver1
"""

import os, sys, os.path
import numpy as np
import pandas as pd
import Processtools as process
from skimage.transform import resize
import Verification as vr

from ftplib import FTP
import netCDF4 as nc


CHIRPSdir = r'E:\CHIRPS-GEFS\Updated_chirps-gefs'

Outdir = r'E:\CHIRPS-GEFS\Process_updated'

DateStart='2016-05-01'
DateEnd = '2019-12-31'
Dates_ =  pd.date_range(DateStart,DateEnd, freq='d')

Dates = Dates_[(Dates_.month >= 6) & (Dates_.month <= 10)]


ftp = FTP('216.218.240.199')
ftp.login(user='ftpuser',passwd=  '@Smekong')
# Files = ftp.nlst('cgefs_precip_a0p05/daily_new/archived/') 

for Date in Dates:
    
    # Date = Dates[0]    
    folder_path = '{0:02d}/{1:02d}/{2:02d}'.format(Date.year,Date.month,Date.day)        
    Files = ftp.nlst('cgefs_precip_a0p05/daily_new/archived/'+ folder_path)     

    Day = 1
    for ftpPath_ in Files:
        # ftpPath_ = Files[0]        
        [Fold,filename] = ftpPath_.split(folder_path+'/')  
        
        MK_files = os.path.join(Outdir,'D{}'.format(Day))      
        chirpsgefs_file = os.path.join(MK_files,filename)     
                
        retry = True
        while (retry):
            try:
            
                with open(chirpsgefs_file, 'wb') as outfile:
                    ftp.retrbinary("RETR " + ftpPath_, outfile.write)       
                    
                retry = False
    
            except IOError as e:
                print( "I/O error({0}): {1}".format(e.errno, e.strerror))
                print( "Retrying...")
                ftp = FTP('216.218.240.199')
                ftp.login(user='ftpuser',passwd=  '@Smekong')
                retry = True  
               
        Day = Day +1 
        
    print(Date)
           


    
    