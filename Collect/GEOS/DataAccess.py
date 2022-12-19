# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Tue Feb 19 08:47:49 2019
"""

import os
import numpy as np
import pandas as pd
import requests
import urllib
import time

def DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, TimeStep, Period, Waitbar, data_type = ["mean"]):

	# watertools modules
    import watertools.General.data_conversions as DC
	
    # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -90 or latlim[1] > 90:
        print('Latitude above 90N or below 90S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -90)
        latlim[1] = np.min(latlim[1], 90)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)  
    
    # Get information of the parameter    
    VarInfo = VariablesInfo(TimeStep)
    Parameter = VarInfo.names[Var]
    unit  = VarInfo.units[Var]
    types  = VarInfo.types[Var]
    
    # Create output folder
    output_folder = os.path.join(Dir, "GEOS", Parameter, TimeStep) 
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Define IDs
    IDx = [np.floor((lonlim[0] + 180)/0.3125), np.ceil((lonlim[1] + 180)/0.3125)]
    IDy = [np.floor((latlim[0] + 90)/0.25), np.ceil((latlim[1] + 90)/0.25)]
    
    # Create output geo transform
    Xstart = -180 + 0.3125 * IDx[0]
    Ystart = -90 + 0.25 * IDy[1]
    geo_out = tuple([Xstart, 0.3125, 0, Ystart, 0, -0.25])
    proj = "WGS84"
    
    Dates = pd.date_range(Startdate, Enddate, freq = "D")
   
    
    # Create Waitbar
    if Waitbar == 1:
        import watertools.Functions.Random.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
    
    for Date in Dates:
        
        # Define the IDz
        if TimeStep == "hourly":
            IDz_start = IDz_end = int(((Date - pd.Timestamp("2017-12-01")).days) * 24) + (Period - 1)
            Hour = int((Period - 1))
            output_name_min = output_folder
            output_name_max = output_folder
            if Var in ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp','t10m', 'v10m', 'u10m', 'v50m', 'u50m', 'ts']:
                
                forecast = 0
                
                if "Max_IDz_start" not in locals():
                    response = requests.get("https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_slv_Nx.info")
                    if response.status_code > 300:
                        raise("Error opening the server https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_slv_Nx.info")
                    
                    content = str(response.content)
                    Max_IDz_start = int(content.split("points")[-2].split("(")[-1]) - 1

                    
                if IDz_start > Max_IDz_start:    
                    
                    Date_start_fc = pd.Timestamp("2017-12-01") + pd.DateOffset(hours = Max_IDz_start)
                    print("forecast is used")
                    url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/tavg1_2d_slv_Nx/tavg1_2d_slv_Nx.%d%d%d_00." %(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)
                    IDz_start = IDz_end = int(((Date - pd.Timestamp(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)).days) * 24) + (Period - 1)
                    forecast = 1
                    
                else:
                    url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_slv_Nx."
                    
            if Var in ['swgdn']:
                
                forecast = 0
                
                if "Max_IDz_start" not in locals():
                    response = requests.get("https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_rad_Nx.info")
                    if response.status_code > 300:
                        raise("Error opening the server https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_rad_Nx.info")
                    
                    content = str(response.content)
                    Max_IDz_start = int(content.split("points")[-2].split("(")[-1]) - 1

                if IDz_start > Max_IDz_start:    
                    
                    Date_start_fc = pd.Timestamp("2017-12-01") + pd.DateOffset(hours = Max_IDz_start)
                    print("forecast is used")
                    url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/tavg1_2d_rad_Nx/tavg1_2d_rad_Nx.%d%d%d_00." %(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)
                    IDz_start = IDz_end = int(((Date - pd.Timestamp(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)).days) * 24) + (Period - 1)
                    forecast = 1
                    
                else:
                    url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_rad_Nx."             
  

            if forecast == 0:
                output_name = os.path.join(output_folder, "%s_GEOS_%s_hourly_%d.%02d.%02d_H%02d.M00.tif"%(Var, unit, Date.year, Date.month, Date.day, Hour))
            else:
                output_name = os.path.join(output_folder, "%s_GEOS-fc_%s_hourly_%d.%02d.%02d_H%02d.M00.tif"%(Var, unit, Date.year, Date.month, Date.day, Hour))
     
        
        if TimeStep == "three_hourly":
            
            IDz_start = IDz_end = int(((Date - pd.Timestamp("2017-12-01")).days) * 8) + (Period - 1)
            Hour = int((Period - 1) * 3)
            forecast = 0
            
            if "Max_IDz_start_3h" not in locals():
               response = requests.get("https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/inst3_2d_asm_Nx.info")
               if response.status_code > 300:
                   raise("Error opening the server https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/inst3_2d_asm_Nx.info")
                    
               content = str(response.content)
               Max_IDz_start_3h = int(content.split("points")[-2].split("(")[-1]) - 1
               
            if IDz_start > Max_IDz_start_3h:    
               
               Date_start_fc = pd.Timestamp("2017-12-01") + pd.DateOffset(hours = Max_IDz_start_3h * 3)
               if Var in ['t2m', 't10m', 'ps', 'qv2m', 'gwettop']:
                   print("forecast is used")
                   url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/inst3_2d_asm_Nx/inst3_2d_asm_Nx.%d%d%d_00." %(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)
                   IDz_start = IDz_end = int(((Date - pd.Timestamp(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)).days) * 24) + (Period - 1)
                   forecast = 1
               else:
                   print("parameter not available in three hourly in forecast mode")
               
            else:
               url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/inst3_2d_asm_Nx."
                       
            
            if forecast == 0:
                output_name = os.path.join(output_folder, "%s_GEOS_%s_3-hourly_%d.%02d.%02d_H%02d.M00.tif"%(Var, unit, Date.year, Date.month, Date.day, Hour))
            else:
                output_name = os.path.join(output_folder, "%s_GEOS-fc_%s_3-hourly_%d.%02d.%02d_H%02d.M00.tif"%(Var, unit, Date.year, Date.month, Date.day, Hour))
                
            output_name_min = output_folder
            output_name_max = output_folder
            
            
        if TimeStep == "daily":
    
            if Var in ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp','t10m', 'v10m', 'u10m', 'v50m', 'u50m', 'ts']:
                IDz_start = int(((Date - pd.Timestamp("2017-12-01")).days) * 24) 
                IDz_end = IDz_start + 23   
                forecast = 0
                if "Max_IDz_start" not in locals():
                    
                    response = requests.get("https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_slv_Nx.info")
                    if response.status_code > 300:
                        raise("Error opening the server https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_slv_Nx.info")
                        
                    content = str(response.content)
                    Max_IDz_start = int(content.split("points")[-2].split("(")[-1]) - 1
                
                if IDz_end > Max_IDz_start:    
                    
                    Date_start_fc = pd.Timestamp("2017-12-01") + pd.DateOffset(hours = Max_IDz_start)
                    print("forecast is used")
                    url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/tavg1_2d_slv_Nx/tavg1_2d_slv_Nx.%d%d%d_00." %(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)
                    IDz_start = int(((Date - pd.Timestamp(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)).days) * 24)
                    IDz_end = IDz_start + 23   
                    forecast = 1
                    
                else:
                    url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_slv_Nx."
                    
                size_z = 24
                
            else:
                
                IDz_start = int(((Date - pd.Timestamp("2017-12-01")).days) * 24) 
                IDz_end = IDz_start + 23  
                forecast = 0
                
                if "Max_IDz_start" not in locals():
                    
                     response = requests.get("https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_rad_Nx.info")
                     if response.status_code > 300:
                         raise("Error opening the server https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_rad_Nx.info")
                         
                     content = str(response.content)
                     Max_IDz_start = int(content.split("points")[-2].split("(")[-1]) - 1
                     
                if IDz_end > Max_IDz_start:    
                     
                     Date_start_fc = pd.Timestamp("2017-12-01") + pd.DateOffset(hours = Max_IDz_start)
                     print("forecast is used")
                     url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/tavg1_2d_rad_Nx/tavg1_2d_rad_Nx.%d%d%d_00." %(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)
                     IDz_start = IDz_end = int(((Date - pd.Timestamp(Date_start_fc.year, Date_start_fc.month, Date_start_fc.day)).days) * 24) + (Period - 1)
                     forecast = 1
                     
                else:
                     url_start = r"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_rad_Nx."             

                size_z = 24
                 
            if "mean" in data_type:
                if forecast == 0:
                    output_name = os.path.join(output_folder, "%s_GEOS_%s_daily_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))
                else:
                    output_name = os.path.join(output_folder, "%s_GEOS-fc_%s_daily_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))
            else:
                output_name = output_folder
                
            if "min" in data_type:
                output_name_min = os.path.join(output_folder, "min", "%smin_GEOS_%s_daily_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))
            else:
                output_name_min = output_folder
            if "max" in data_type:
                output_name_max = os.path.join(output_folder, "max", "%smax_GEOS_%s_daily_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))
            else:
                output_name_max = output_folder
                
        if not (os.path.exists(output_name) and os.path.exists(output_name_min) and os.path.exists(output_name_max)):

            # define total url
            url_GEOS = url_start + 'ascii?%s[%s:1:%s][%s:1:%s][%s:1:%s]' %(Var, IDz_start,IDz_end, int(IDy[0]),int(IDy[1]),int(IDx[0]),int(IDx[1]))
        
            # Reset the begin parameters for downloading
            downloaded = 0
            N = 0
            
            # if not downloaded try to download file
            while downloaded == 0:
                try:
    
                    # download data (first save as text file)
                    pathtext = os.path.join(output_folder,'temp%s.txt' %str(IDz_start))
                    
                    # Download the data
                    #print(url_GEOS)
                    urllib.request.urlretrieve(url_GEOS, filename=pathtext)
    
                    # Reshape data
                    datashape = [int(IDy[1] - IDy[0] + 1), int(IDx[1] - IDx[0] + 1)]
                    data_start = np.genfromtxt(pathtext,dtype = float,skip_header = 1,skip_footer = 6,delimiter = ',')
                    data_list = np.asarray(data_start[:,1:])
                    if TimeStep == "daily":
                        data_end = np.resize(data_list,(size_z, datashape[0], datashape[1]))
                    if TimeStep == "hourly" or TimeStep == "three_hourly":
                        data_end = np.resize(data_list,(datashape[0], datashape[1]))                        
                    os.remove(pathtext)
    
                    # Set no data value
                    data_end[data_end>1000000] = np.nan
                    
                    if TimeStep == "daily":

                        if "min" in data_type:
                            data_min = np.nanmin(data_end, 0)                            
                        if "max" in data_type:    
                            data_max = np.nanmax(data_end, 0)                        
                        if "mean" in data_type:
                            if types == "state":
                                data_end = np.nanmean(data_end, 0)
                            else:
                                data_end = np.nansum(data_end, 0)                        
    
                    # Add the VarFactor
                    if VarInfo.factors[Var] < 0:
                        if "mean" in data_type:                        
                            data_end[data_end != -9999] = data_end[data_end != -9999] + VarInfo.factors[Var]
                            data_end[data_end < -9999] = -9999
                            data_end = np.flipud(data_end)
                        if "min" in data_type:
                            data_min[data_min != -9999] = data_min[data_min != -9999] + VarInfo.factors[Var]
                            data_min[data_min < -9999] = -9999
                            data_min = np.flipud(data_min)                                 
                        if "max" in data_type:
                            data_max[data_max != -9999] = data_max[data_max != -9999] + VarInfo.factors[Var]      
                            data_max[data_max < -9999] = -9999
                            data_max = np.flipud(data_max)                        
                                                    
                    else:
                        if "mean" in data_type: 
                            data_end[data_end != -9999] = data_end[data_end != -9999] * VarInfo.factors[Var]
                            data_end[data_end < -9999] = -9999
                            data_end = np.flipud(data_end)
                        if "min" in data_type:
                            data_min[data_min != -9999] = data_min[data_min != -9999] * VarInfo.factors[Var]
                            data_min[data_min < -9999] = -9999
                            data_min = np.flipud(data_min)                            
                        if "max" in data_type:
                            data_max[data_max != -9999] = data_max[data_max != -9999] * VarInfo.factors[Var]                          
                            data_max[data_max < -9999] = -9999
                            data_max = np.flipud(data_max)   
    
                    if "mean" in data_type: 
                        DC.Save_as_tiff(output_name, data_end, geo_out, proj)
                    if "min" in data_type:
                        DC.Save_as_tiff(output_name_min, data_min, geo_out, proj)
                    if "max" in data_type:
                        DC.Save_as_tiff(output_name_max, data_max, geo_out, proj)
                        
                    # Download was succesfull
                    downloaded = 1                        
            
                # If download was not succesfull
                except Exception as e:
                    print(e)
                    # Try another time
                    N = N + 1
    
                    # Stop trying after 10 times
                    if N == 10:
                        print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                        time.sleep(30)
                        downloaded = 1
    
        if Waitbar == 1:
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    return()    

class VariablesInfo:
    """
    This class contains the information about the GEOS variables
    """
    names = {'t2m': 'Air_Temperature',
             'u2m': 'Eastward_Wind',
             'v2m': 'Northward_Wind',
             'qv2m': 'Specific_Humidity',
             'tqv': 'Total_Precipitable_Water_Vapor',
             'ps': 'Surface_Pressure',
             'slp': 'Sea_Level_Pressure',
             't10m': 'Air_Temperature_10m',
             'v10m': 'Northward_Wind_10m',
             'u10m': 'Eastward_Wind_10m',
             'v50m': 'Northward_Wind_50m',
             'u50m': 'Eastward_Wind_50m',
             'ts': 'Surface_Skin_Temperature',
             'emis': 'Surface_Emissivity',
             'swgnt': 'Surface_Net_Downward_Shortwave_Flux',
             'swgntclr': 'Surface_Net_Downward_Shortwave_Flux_Assuming_Clear_Sky',
             'swgntcln': 'Surface_Net_Downward_Shortwave_Flux_Assuming_No_Aerosol',         
             'swgntclrcln': 'Surface_Net_Downward_Shortwave_Flux_Assuming_Clear_Sky_And_No_Aerosol',
             'swgdn': 'Surface_Incoming_Shortwave_Flux',
             'lwgnt': 'Surface_Net_Downward_Longwave_Flux',
             'lwgntclr': 'Surface_Net_Downward_Longwave_Flux_Assuming_Clear_Sky',
             'lwgntclrcln': 'Surface_Net_Downward_Longwave_Flux_Assuming_Clear_Sky_And_No_Aerosol'}
    
    descriptions = {'t2m': '2m Air Temperature',
             'u2m': '2m Eastward wind',
             'v2m': '2m Northward wind',
             'qv2m': '2m Specific Humidity',
             'tqv': 'Total Precipitable Water Vapor',
             'ps': 'Surface Pressure',
             'slp': 'Sea Level Pressure',
             't10m': '10m Air Temperature',
             'v10m': '10m Northward wind',
             'u10m': '10m Eastward wind',
             'v50m': '50m Northward wind',
             'u50m': '50m Eastward wind',
             'ts': 'Surface Skin Temperature',
             'emis': 'Surface Emissivity',
             'swgnt': 'Surface Net Downward Shortwave Flux',
             'swgntclr': 'Surface Net Downward Shortwave Flux Assuming Clear Sky',
             'swgntcln': 'Surface Net Downward Shortwave Flux Assuming No Aerosol',         
             'swgntclrcln': 'Surface Net Downward Shortwave Flux Assuming Clear Sky And No Aerosol',
             'swgdn': 'Surface Incoming Shortwave Flux',
             'lwgnt': 'Surface Net Downward Longwave Flux',
             'lwgntclr': 'Surface Net Downward Longwave Flux Assuming Clear Sky',
             'lwgntclrcln': 'Surface Net Downward Longwave Flux Assuming Clear Sky And No Aerosol'}
    
    factors = {'t2m': 1,
             'u2m': 1,
             'v2m': 1,
             'qv2m': 1,
             'tqv': 1,
             'ps': 0.001,
             'slp':  0.001,
             't10m': 1,
             'v10m': 1,
             'u10m': 1,
             'v50m': 1,
             'u50m': 1,
             'ts': 1,
             'emis': 1,
             'swgnt': 1,
             'swgntclr': 1,
             'swgntcln': 1,         
             'swgntclrcln': 1,
             'swgdn': 1,
             'lwgnt': 1,
             'lwgntclr': 1,
             'lwgntclrcln': 1}
    
    types = {'t2m': 'state',
             'u2m': 'state',
             'v2m': 'state',
             'qv2m': 'state',
             'tqv': 'state',
             'ps': 'state',
             'slp': 'state',
             't10m': 'state',
             'v10m': 'state',
             'u10m': 'state',
             'v50m': 'state',
             'u50m': 'state',
             'ts': 'state',
             'emis': 'state',
             'swgnt': 'state',
             'swgntclr': 'state',
             'swgntcln': 'state',         
             'swgntclrcln': 'state',
             'swgdn': 'state',
             'lwgnt': 'state',
             'lwgntclr': 'state',
             'lwgntclrcln': 'state'}

    def __init__(self, step):
        if step == 'three_hourly':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'qv2m': 'kg-kg-1',
             'tqv': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             't10m': 'K',
             'v10m': 'm-s-1',
             'u10m': 'm-s-1',
             'v50m': 'm-s-1',
             'u50m': 'm-s-1'}
            
        elif step == 'hourly':
            self.units = {'ts': 'K',
             'emis': '-',
             'swgnt': 'W-m-2',
             'swgntclr': 'W-m-2',
             'swgntcln': 'W-m-2',         
             'swgntclrcln': 'W-m-2',
             'swgdn': 'W-m-2',
             'lwgnt': 'W-m-2',
             'lwgntclr': 'W-m-2',
             'lwgntclrcln': 'W-m-2',
             't2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'qv2m': 'kg-kg-1',
             'tqv': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             't10m': 'K',
             'v10m': 'm-s-1',
             'u10m': 'm-s-1',
             'v50m': 'm-s-1',
             'u50m': 'm-s-1'}  
            
            
        elif step == 'daily':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'qv2m': 'kg-kg-1',
             'tqv': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             't10m': 'K',
             'v10m': 'm-s-1',
             'u10m': 'm-s-1',
             'v50m': 'm-s-1',
             'u50m': 'm-s-1',
             'ts': 'K',
             'emis': '-',
             'swgnt': 'W-m-2',
             'swgntclr': 'W-m-2',
             'swgntcln': 'W-m-2',         
             'swgntclrcln': 'W-m-2',
             'swgdn': 'W-m-2',
             'lwgnt': 'W-m-2',
             'lwgntclr': 'W-m-2',
             'lwgntclrcln': 'W-m-2'}

        else:
            raise KeyError("The input time step is not supported")

