# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:24:54 2020

@author: timhe
"""

# General Python modules
import numpy as np
import os
import glob
import pandas as pd
import gdal

def Nearest_Interpolate(Dir_in, Startdate, Enddate, Dir_out=None):
    """
    This functions calculates monthly tiff files based on the daily tiff files.
    (will calculate the total sum)

    Parameters
    ----------
    Dir_in : str
        Path to the input data
    Startdate : str
        Contains the start date of the model 'yyyy-mm-dd'
    Enddate : str
        Contains the end date of the model 'yyyy-mm-dd'
    Dir_out : str
        Path to the output data, default is same as Dir_in

    """
    # import WA+ modules
    import watertools.General.data_conversions as DC
    import watertools.General.raster_conversions as RC

    # Change working directory
    os.chdir(Dir_in)

    # Define end and start date
    Dates = pd.date_range(Startdate, Enddate, freq='D')

    # Find all hourly files
    files = glob.glob('*hourly*.tif')

    # Get array information and define projection
    geo_out, proj, size_X, size_Y = RC.Open_array_info(files[0])
    if int(proj.split('"')[-2]) == 4326:
        proj = "WGS84"

    # Get the No Data Value
    dest = gdal.Open(files[0])
    NDV = dest.GetRasterBand(1).GetNoDataValue()

    # Define output directory
    if Dir_out is None:
	     Dir_out = Dir_in

    if not os.path.exists(Dir_out):
	     os.makedirs(Dir_out)

    # loop over the months and sum the days
    for date in Dates:
        Year = date.year
        Month = date.month
        Day = date.day
        files_one_day = glob.glob('*hourly*%d.%02d.%02d*.tif' % (Year, Month, Day))

        # Create empty arrays
        Daily_data = np.zeros([size_Y, size_X])

        if len(files_one_day) != 24:
            print("One hour is missing!!! day %s month %s year %s" %(Day, Month, Year))

        for file_one_year in files_one_day:
            file_path = os.path.join(Dir_in, file_one_year)

            Hour_data = RC.Open_tiff_array(file_path)
            Hour_data[np.isnan(Hour_data)] = 0.0
            Hour_data[Hour_data == -9999] = 0.0
            Daily_data += Hour_data

        # Define output name
        output_name = os.path.join(Dir_out, file_one_year
                                   .replace('hourly', 'daily')
                                   .replace('hour', 'day'))

        output_name = output_name[:-19] + '%d.%02d.%02d.tif' % (date.year, date.month, date.day)

        # Save tiff file
        DC.Save_as_tiff(output_name, Daily_data, geo_out, proj)

    return
