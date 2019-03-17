# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Sun Feb 10 18:26:30 2019
"""
import datetime
import pandas as pd
import numpy as np
import os
import urllib
import gdal

#S2_tile = "10SGE"
#output_folder = "F:\Project_Jain\Case_California\Input_Data\S30"
#Startdate = "2018-03-01"
#Enddate = "2018-10-31"

def DownloadData(Dir, Startdate, Enddate, S2_tile):

    # Import watertools modules
    import watertools.General.data_conversions as DC
    
    # Define the dates
    Dates = pd.date_range(Startdate, Enddate, freq = "D")
    
    # Create output folder
    output_folder_end = os.path.join(Dir, S2_tile)
    if not os.path.exists(output_folder_end):
        os.makedirs(output_folder_end)
    
    # Loop over the days
    for Date in Dates:
        
        # Get the datum
        doy = int(Date.dayofyear)
        year = Date.year
        
        try:
            
            # Create the right url
            url = "https://hls.gsfc.nasa.gov/data/v1.4/S30/%s/%s/%s/%s/%s/HLS.S30.T%s.%s%03d.v1.4.hdf" % (year, S2_tile[0:2], S2_tile[2:3], S2_tile[3:4], S2_tile[4:5],S2_tile, year, doy)
            filename_out = os.path.join(output_folder_end, "HLS.S30.T%s.%s%03d.v1.4.hdf" % (S2_tile, year, doy))

            # Download the data
            urllib.request.urlretrieve(url, filename=filename_out)
            
            # Create a folder for the end results 
            folder_tile = os.path.join(output_folder_end, "HLS_S30_T%s_%s%03d_v1_4" % (S2_tile, year, doy))
            if not os.path.exists(folder_tile):
                os.makedirs(folder_tile)        
            
            # Write the hdf file in tiff files and store it in the folder
            dataset = gdal.Open(filename_out)
            sdsdict = dataset.GetMetadata('SUBDATASETS')
            
            for Band in range(1,15):
                dest = gdal.Open(sdsdict["SUBDATASET_%d_NAME" %Band])
                Array = dest.GetRasterBand(1).ReadAsArray()
                if Band < 9.:
                    Array = Array * 0.0001
                    Array[Array == -0.1] = -9999.
                    Band_name = "B%02d" %(Band)
                if Band == 9.:
                    Band_name = "B8A"
                    Array = Array * 0.0001
                    Array[Array == -0.1] = -9999.
                if (Band >= 10. and Band < 14.):
                    Band_name = "B%02d" %(Band - 1)   
                    Array = Array * 0.0001
                    Array[Array == -0.1] = -9999.
                if Band == 14.:
                    Array[Array == 255] = -9999.
                    Band_name = "QC"
                
                meta = dataset.GetMetadata()
                ulx = int(meta["ULX"])
                uly = int(meta["ULY"])   
                size = int(meta["SPATIAL_RESOLUTION"]) 
                projection = int(meta["HORIZONTAL_CS_CODE"].split(":")[-1])
                
                time_string = meta["SENSING_TIME"].split(";")[0]

                Time = datetime.datetime.strptime(time_string[:-5],"%Y-%m-%dT%H:%M:%S")
                hour = int(Time.hour)
                minute = int(Time.minute)
                geo = tuple([ulx, size, 0, uly, 0, -size])
                DC.Save_as_tiff(os.path.join(folder_tile, "HLS_S30_T%s_%s%03d_H%02dM%02d_%s.tif" % (S2_tile, year, doy, hour, minute, Band_name)), Array, geo, projection)
                     
        except:
            pass
        
    return()