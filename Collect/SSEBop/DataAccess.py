# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/SSEBop

Restrictions:
The data and this python file may not be distributed to others without
permission of the WA+ team due data restriction of the SSEBop developers.

Description:
This script collects SSEBop data from the UNESCO-IHE FTP server or from the web. The data has a
monthly temporal resolution and a spatial resolution of 0.01 degree. The
resulting tiff files are in the WGS84 projection.
The data is available between 2003-01-01 till present.

Example:
from watertools.Collect import SSEBop
SSEBop.monthly(Dir='C:/Temp/', Startdate='2003-02-24', Enddate='2003-03-09',
                     latlim=[50,54], lonlim=[3,7])

"""
# General modules
import numpy as np
import os
import pandas as pd
from ftplib import FTP
import sys
if sys.version_info[0] == 3:
    import urllib.parse
if sys.version_info[0] == 2:
    import urllib
# Water Accounting Modules
import watertools
import watertools.General.raster_conversions as RC
import watertools.General.data_conversions as DC


def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar, version, TimeStep, Product):
    """
    This scripts downloads SSEBop ET data from the UNESCO-IHE ftp server.
    The output files display the total ET in mm for a period of one month.
    The name of the file corresponds to the first day of the month.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    lonlim -- [ymin, ymax] (values must be between -90 and 90)
    latlim -- [xmin, xmax] (values must be between -180 and 180)
    """

    if version == "FTP":
        # Check the latitude and longitude and otherwise set lat or lon on greatest extent
        if latlim[0] < -59.2 or latlim[1] > 80:
            print('Latitude above 80N or below -59.2S is not possible. Value set to maximum')
            latlim[0] = np.max(latlim[0], -59.2)
            latlim[1] = np.min(latlim[1], 80)
        if lonlim[0] < -180 or lonlim[1] > 180:
            print('Longitude must be between 180E and 180W. Now value is set to maximum')
            lonlim[0] = np.max(lonlim[0],-180)
            lonlim[1] = np.min(lonlim[1],180)

    	# Check Startdate and Enddate
        if not Startdate:
            Startdate = pd.Timestamp('2003-01-01')
        if not Enddate:
            Enddate = pd.Timestamp('2014-10-31')

    if version == "V4":
        # Check the latitude and longitude and otherwise set lat or lon on greatest extent
        if latlim[0] < -60 or latlim[1] > 80.0022588483988670:
            print('Latitude above 80N or below -59.2S is not possible. Value set to maximum')
            latlim[0] = np.max(latlim[0], -60)
            latlim[1] = np.min(latlim[1], 80.0022588483988670)
        if lonlim[0] < -180 or lonlim[1] > 180.0002930387853439:
            print('Longitude must be between 180E and 180W. Now value is set to maximum')
            lonlim[0] = np.max(lonlim[0],-180)
            lonlim[1] = np.min(lonlim[1],180.0002930387853439)

    	# Check Startdate and Enddate
        if not Startdate:
            Startdate = pd.Timestamp('2003-01-01')
        if not Enddate:
            import datetime
            Enddate = pd.Timestamp(datetime.datetime.now())

   # Define directory and create it if not exists
    if TimeStep == "daily":
        if Product == "ETact":
            output_folder = os.path.join(Dir, 'Evapotranspiration', 'SSEBop', 'Daily')
            freq_use = "D"
        if Product == "ETpot":
            output_folder = os.path.join(Dir, 'Potential_Evapotranspiration', 'FEWS', 'Daily')        
            freq_use = "D"

    if TimeStep == "monthly":
        if Product == "ETact":
            output_folder = os.path.join(Dir, 'Evapotranspiration', 'SSEBop', 'Monthly')
            freq_use = "MS"
        if Product == "ETpot":
            output_folder = os.path.join(Dir, 'Potential_Evapotranspiration', 'FEWS', 'Monthly')        
            freq_use = "MS"

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Creates dates library
    Dates = pd.date_range(Startdate, Enddate, freq = freq_use)

    # Create Waitbar
    if Waitbar == 1:
        import watertools.Functions.Random.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    # Loop over the dates
    for Date in Dates:

        # Define year and month
        year = Date.year
        month = Date.month
        day = Date.day

        if version == "V4" or version == "V5" or version == "V6":

            # Date as printed in filename
            if Product == "ETpot":
                
                if TimeStep == "daily":
                    Filename_out= os.path.join(output_folder,'ETpot_FEWS_mm-day-1_daily_%s.%02s.%02s.tif' %(Date.strftime('%Y'), Date.strftime('%m'), Date.strftime('%d')))
                    # Define the downloaded zip file
                    Filename_only_zip = 'et%02s%02d%02d.tar.gz' %(str(year)[2:], month, day)
                    # The end file name after downloading and unzipping
                    Filename_only = "et%02s%02d%02d.bil" %(str(year)[2:], month, day)
                    # Create bin folder
                    temp_folder = os.path.join(output_folder, "Temp")
                    if not os.path.exists(temp_folder):
                        os.makedirs(temp_folder)
                    local_filename = os.path.join(temp_folder, Filename_only)
                    
            if Product == "ETact" and version == "V4":
                                    
                if TimeStep == "monthly":
                    Filename_out= os.path.join(output_folder,'ETa_SSEBop_V4_mm-month-1_monthly_%s.%02s.%02s.tif' %(Date.strftime('%Y'), Date.strftime('%m'), Date.strftime('%d')))
                    # Define the downloaded zip file
                    Filename_only_zip = "m%s%02d.zip" %(str(year), month)
                    # The end file name after downloading and unzipping
                    Filename_only = "m%s%02d_modisSSEBopETv4_actual_mm.tif" %(str(year), month)
    
            		    # Temporary filename for the downloaded global file
                    local_filename = os.path.join(output_folder, Filename_only)

            if Product == "ETact" and version == "V5":
                
                if TimeStep == "monthly":
                    Filename_out= os.path.join(output_folder,'ETa_SSEBop_V5_mm-month-1_monthly_%s.%02s.%02s.tif' %(Date.strftime('%Y'), Date.strftime('%m'), Date.strftime('%d')))
                    # Define the downloaded zip file
                    Filename_only_zip = "m%s%02d.zip" %(str(year), month)
                    # The end file name after downloading and unzipping
                    Filename_only = "m%s%02d_modisSSEBopETv5_actual_mm.tif" %(str(year), month)
    
            		    # Temporary filename for the downloaded global file
                    local_filename = os.path.join(output_folder, Filename_only)

            if Product == "ETact" and version == "V6":
                
                if TimeStep == "monthly":
                    Filename_out= os.path.join(output_folder,'ETa_SSEBop_V6_mm-month-1_monthly_%s.%02s.%02s.tif' %(Date.strftime('%Y'), Date.strftime('%m'), Date.strftime('%d')))
                    # Define the downloaded zip file
                    Filename_only_zip = "m%s%02d.zip" %(str(year), month)
                    # The end file name after downloading and unzipping
                    Filename_only = "m%s%02d_viirsSSEBopETv6_actual_mm.tif" %(str(year), month)
    
            		    # Temporary filename for the downloaded global file
                    local_filename = os.path.join(output_folder, Filename_only)




        # Download the data from FTP server if the file not exists
        if not os.path.exists(Filename_out):
            try:

                if version == "V4" or version == "V5" or version == "V6":
                    if Product == "ETpot":
                        Download_SSEBop_from_Web(temp_folder, Filename_only_zip, Product, TimeStep, version)
                    if Product == "ETact":
                        Download_SSEBop_from_Web(output_folder, Filename_only_zip, Product, TimeStep, version)
                        
                if Product == "ETpot":
                    Array_ETpot = RC.Open_bil_array(local_filename)
                    Array_ETpot = Array_ETpot/100
                    Geo_out = tuple([-180.5, 1, 0, 90.5, 0, -1])
                    dest = DC.Save_as_MEM(Array_ETpot, Geo_out, "WGS84")
                    data, Geo_out, proj = RC.clip_data(dest, latlim, lonlim)
                    DC.Save_as_tiff(Filename_out, data, Geo_out, "WGS84")
                    
                if Product == "ETact":    
                    # Clip dataset
                    data, Geo_out, proj = RC.clip_data(local_filename, latlim, lonlim)
                    data[data<-9999] = -9999
                    DC.Save_as_tiff(Filename_out, data, Geo_out, "WGS84")
                    os.remove(local_filename)

            except:
                print("Was not able to download file with date %s" %Date)

        # Adjust waitbar
        if Waitbar == 1:
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    if version == "V4" or version == "V5":
        import glob
        os.chdir(output_folder)
        if Product == "ETact":
            zipfiles = glob.glob("*.zip")
            for zipfile in zipfiles:
                os.remove(os.path.join(output_folder, zipfile))
            xmlfiles = glob.glob("*.xml")
            for xmlfile in xmlfiles:
                os.remove(os.path.join(output_folder, xmlfile))
        if Product == "ETpot":  
            import shutil
            Temp_dir = os.path.join(output_folder, "Temp")
            shutil.rmtree(Temp_dir)
            
    return

def Download_SSEBop_from_Web(output_folder, Filename_only_zip, Product, TimeStep, version):
    """
    This function retrieves SSEBop data for a given date from the
    https://edcintl.cr.usgs.gov server.

    Keyword arguments:
	 local_filename -- name of the temporary file which contains global SSEBop data
    Filename_dir -- name of the end directory to put in the extracted data
    """
    if Product == "ETact" and TimeStep == "monthly" and version == "V4":
        # Create the total url to the webpage
        total_URL = "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/fews/web/global/monthly/eta/downloads/" + str(Filename_only_zip)

    if Product == "ETact" and TimeStep == "monthly" and version == "V5":
        # Create the total url to the webpage
        total_URL = "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/fews/web/global/monthly/etav5/downloads/" + str(Filename_only_zip)

    if Product == "ETact" and TimeStep == "monthly" and version == "V6":
        # Create the total url to the webpage
        total_URL = "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/fews/web/global/monthly/etav6/downloads/monthly/" + str(Filename_only_zip)

           
    if Product == "ETpot" and TimeStep == "daily":
        total_URL = "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/fews/web/global/daily/pet/downloads/daily/" + str(Filename_only_zip)



    # Download the data
    if sys.version_info[0] == 2:
        urllib.urlretrieve(total_URL, os.path.join(output_folder, Filename_only_zip))
    if sys.version_info[0] == 3:
        urllib.request.urlretrieve(total_URL, os.path.join(output_folder, Filename_only_zip)) 

    # unzip the file
    if Product == "ETpot":
        DC.Extract_Data_tar_gz(os.path.join(output_folder, Filename_only_zip), output_folder)
    if Product == "ETact":
        DC.Extract_Data(os.path.join(output_folder, Filename_only_zip), output_folder)

    return

