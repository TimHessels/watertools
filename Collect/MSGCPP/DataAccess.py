# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Tue Feb 19 07:10:44 2019
"""
import os
import requests
import datetime
import pandas as pd
import numpy as np
from netCDF4 import Dataset

import watertools.General.data_conversions as DC

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Time = '', GMT_Offset = 0, Waitbar = 1, Type = "SDS"):
    
    # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -90 or latlim[1] > 90:
        print('Latitude above 90N or below 90S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -90)
        latlim[1] = np.min(latlim[1], 90)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)
    
    
    output_folder = os.path.join(Dir, "MSGCPP", Type, "15min")
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
            
        if Type == 'SDS':
            filename_out = os.path.join(output_folder, "SDS_MSGCPP_W-m-2_15min_%d.%02d.%02d_H%02d.M%02d.tif" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))
            filename_out_nc = os.path.join(output_folder, "SDS_MSGCPP_W-m-2_15min_%d.%02d.%02d_H%02d.M%02d.nc" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))
            parameter = 'sds'

        if Type == 'lwe_precipitation_rate':
            filename_out = os.path.join(output_folder, "Precipitation_MSGCPP_mm-h-1_15min_%d.%02d.%02d_H%02d.M%02d.tif" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))
            filename_out_nc = os.path.join(output_folder, "Precipitation_MSGCPP_mm-h-1_15min_%d.%02d.%02d_H%02d.M%02d.nc" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))
            parameter = 'precip'
            
        if Type == 'lwe_precipitation_rate_ir':
            filename_out = os.path.join(output_folder, "PrecipitationIR_MSGCPP_mm-h-1_15min_%d.%02d.%02d_H%02d.M%02d.tif" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))
            filename_out_nc = os.path.join(output_folder, "PrecipitationIR_MSGCPP_mm-h-1_15min_%d.%02d.%02d_H%02d.M%02d.nc" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))
            parameter = 'precip_ir'
            
        if Type == 'Cloud':
            filename_out = os.path.join(output_folder, "Cloud_MSGCPP_-_15min_%d.%02d.%02d_H%02d.M%02d.tif" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))
            filename_out_nc = os.path.join(output_folder, "Cloud_MSGCPP_-_15min_%d.%02d.%02d_H%02d.M%02d.nc" %(Date.year, Date.month, Date.day, Date.hour, Date.minute))
            parameter = 'cldmask'

        if not os.path.exists(filename_out):
        
            # define url
            if Type == 'SDS':
                url = r"https://msgcpp-adaguc.knmi.nl/adaguc-server?dataset=msgrt&service=wcs&request=getcoverage&coverage=surface_downwelling_shortwave_flux_in_air&FORMAT=NetCDF4&CRS=EPSG%%3A4326&BBOX=%s,%s,%s,%s&RESX=0.04310344827586207&RESY=0.04418103448275862&time=%d-%02d-%02dT%02d%%3A%02d%%3A00Z" %(lonlim[0],latlim[0], lonlim[1], latlim[1], Date.year, Date.month, Date.day, Date.hour, Date.minute)
            if Type == 'lwe_precipitation_rate':
                url = r"https://msgcpp-adaguc.knmi.nl/adaguc-server?dataset=msgrt&service=wcs&request=getcoverage&coverage=lwe_precipitation_rate&FORMAT=NetCDF4&CRS=EPSG%%3A4326&BBOX=%s,%s,%s,%s&RESX=0.04310344827586207&RESY=0.04418103448275862&time=%d-%02d-%02dT%02d%%3A%02d%%3A00Z" %(lonlim[0],latlim[0], lonlim[1], latlim[1], Date.year, Date.month, Date.day, Date.hour, Date.minute)
            if Type == 'lwe_precipitation_rate_ir':
                url = r"https://msgcpp-adaguc.knmi.nl/adaguc-server?dataset=msgrt&service=wcs&request=getcoverage&coverage=lwe_precipitation_rate_ir&FORMAT=NetCDF4&CRS=EPSG%%3A4326&BBOX=%s,%s,%s,%s&RESX=0.04310344827586207&RESY=0.04418103448275862&time=%d-%02d-%02dT%02d%%3A%02d%%3A00Z" %(lonlim[0],latlim[0], lonlim[1], latlim[1], Date.year, Date.month, Date.day, Date.hour, Date.minute)
            if Type == 'Cloud':
                url = r"https://msgcpp-adaguc.knmi.nl/adaguc-server?dataset=msgrt&service=wcs&request=getcoverage&coverage=cloud_area_fraction&FORMAT=NetCDF4&CRS=EPSG%%3A4326&BBOX=%s,%s,%s,%s&RESX=0.04310344827586207&RESY=0.04418103448275862&time=%d-%02d-%02dT%02d%%3A%02d%%3A00Z" %(lonlim[0],latlim[0], lonlim[1], latlim[1], Date.year, Date.month, Date.day, Date.hour, Date.minute)

            print(url)
            success = 0
            attempts = 0
            while success == 0 and attempts<4:
                try:
                    request = requests.get(url, timeout=10, stream=True)
                    with open(filename_out_nc, 'wb') as fh:
                        # Walk through the request response in chunks of 1024 * 1024 bytes, so 1MiB
                        for chunk in request.iter_content(1024 * 1024):
                            # Write the chunk to the file
                            fh.write(chunk)
                            # Optionally we can check here if the download is taking too long                    
                    #urllib.request.urlretrieve(url, filename=)
                    statinfo = os.stat(filename_out_nc)
                    if statinfo.st_size > 300:
                        success = 1
                        
                    else:
                        attempts +=1
                        
                        
                except Exception as e:
                    print("ERROR: %s" %e)     
                    attempts +=1

            if statinfo.st_size < 300:
                os.remove(filename_out_nc)
                
                # maak hier nep file met nans
                #0.04310344827586207&RESY=0.04418103448275862
            else:
                
                fh = Dataset(filename_out_nc)
                Array = np.array(fh[parameter][:,:].data.squeeze())
                Array[Array<0] = np.nan
                x = fh["x"][:]               
                y = fh["y"][:]            
                x_size = (np.nanmax(x) - np.nanmin(x))/(len(x)-1)
                y_size = (np.nanmax(y) - np.nanmin(y))/(len(y)-1)
                lon_start = np.nanmin(x) - x_size/2
                lat_start = np.nanmax(y) + y_size/2
                proj = 4326
                geo = tuple([lon_start, x_size, 0, lat_start, 0, -y_size])
                
                DC.Save_as_tiff(filename_out, Array, geo, proj)
                
                try:
                    fh.close()
                    os.remove(filename_out_nc)
                except Exception as e:
                    print("ERROR removing %s" %filename_out_nc)
                    print("MESSAGE: %s" %e)               
        
        if Waitbar == 1:
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
        