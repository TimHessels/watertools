# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Tue Feb 19 07:10:44 2019
"""
import os
import urllib
import datetime
import pandas as pd
import numpy as np

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Time = '', GMT_Offset = 0, Waitbar = 1):
    
    # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -90 or latlim[1] > 90:
        print('Latitude above 90N or below 90S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -90)
        latlim[1] = np.min(latlim[1], 90)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)
    
    output_folder = os.path.join(Dir, "MSGCPP", "SDS", "15min")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    if isinstance(Enddate, str):
        Enddate = datetime.datetime(int(Enddate.split('-')[0]), int(Enddate.split('-')[1]), int(Enddate.split('-')[2]), 23, 59)
    else:    
        Enddate = datetime.datetime(Enddate.year, Enddate.month, Enddate.day, 23, 59)
    
    if Time == '':
        Dates = pd.date_range(Startdate, Enddate, freq = "15min") - datetime.timedelta(hours=GMT_Offset)
    else:
        Dates = pd.date_range(Startdate, Enddate, freq = "D") - datetime.timedelta(hours=GMT_Offset)
    
    # Create Waitbar
    if Waitbar == 1:
        import watertools.Functions.Random.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
    # Loop over dates
    for Date in Dates:
        
        if Time != '':
            Hour = int(Time.split(':')[0])
            Minute = int(Time.split(':')[1])        
            Date = datetime.datetime(Date.year, Date.month, Date.day, Hour, Minute)
        
        filename_out = os.path.join(output_folder, "SDS_MSGCPP_W-m-2_15min_%d.%02d.%02d_H%02d.M%02d.tif" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))

        if not os.path.exists(filename_out):
        
            # define url
            url = r"http://msgcpp-ogc-archive.knmi.nl/msgar.cgi?&service=wcs&version=1.0.0&request=getcoverage&coverage=surface_downwelling_shortwave_flux_in_air&FORMAT=GeoTIFF&CRS=EPSG%%3A4326&BBOX=%s,%s,%s,%s&RESX=0.04310344827586207&RESY=0.04418103448275862&time=%d-%02d-%02dT%02d%%3A%02d%%3A00Z" %(lonlim[0],latlim[0], lonlim[1], latlim[1], Date.year, Date.month, Date.day, Date.hour, Date.minute)
            urllib.request.urlretrieve(url, filename=filename_out)
            statinfo = os.stat(filename_out)

            if statinfo.st_size < 3000:

                #url = r"https://msgcpp-adaguc.knmi.nl/adaguc-server/?dataset=msgrt&service=wcs&request=getcoverage&coverage=surface_downwelling_shortwave_flux_in_air&FORMAT=GeoTIFF&CRS=EPSG%%3A4326&&BBOX=%s,%s,%s,%s&RESX=0.04310344827586207&RESY=0.04418103448275862&time=%d-%02d-%02dT%02d%%3A%02d%%3A00Z"  %(lonlim[0],latlim[0], lonlim[1], latlim[1], Date.year, Date.month, Date.day, Date.hour, Date.minute)
                #print("new URL %s" %url)
                url = r"http://msgcpp-ogc-realtime.knmi.nl/msgrt.cgi??&service=wcs&version=1.0.0&request=getcoverage&coverage=surface_downwelling_shortwave_flux_in_air&FORMAT=GeoTIFF&CRS=EPSG%%3A4326&BBOX=%s,%s,%s,%s&RESX=0.04310344827586207&RESY=0.04418103448275862&time=%d-%02d-%02dT%02d%%3A%02d%%3A00Z" %(lonlim[0],latlim[0], lonlim[1], latlim[1], Date.year, Date.month, Date.day, Date.hour, Date.minute)
                urllib.request.urlretrieve(url, filename=filename_out)         
                statinfo = os.stat(filename_out)

            if statinfo.st_size < 300:
                os.remove(filename_out)
                
                # maak hier nep file met nans
                #0.04310344827586207&RESY=0.04418103448275862
           
            
        
        if Waitbar == 1:
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
        