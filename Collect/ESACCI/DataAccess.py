# -*- coding: utf-8 -*-

# General modules
import shutil
import os
import urllib

# Water Accounting modules
import watertools.General.data_conversions as DC
import watertools.General.raster_conversions as RC

def DownloadData(Dir, latlim, lonlim, Waitbar):
    """
    This function downloads NLDAS Forcing data hourly, daily or monthly data

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    """
    
    # Define the output name
    output_filename = os.path.join(Dir, 'LU_ESACCI.tif')
    
    # Set the url of the server
    url = r"https://storage.googleapis.com/cci-lc-v207/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.zip"
    
    # Create a Trash folder
    Dir_trash = os.path.join(Dir, "Trash")
    if not os.path.exists(Dir_trash):
        os.makedirs(Dir_trash)
    
    # Define location of download
    filename_out = os.path.join(Dir_trash, "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.zip")

    # Download data
    urllib.request.urlretrieve(url, filename=filename_out)
    
    # Extract data
    DC.Extract_Data(filename_out, Dir_trash)
    
    # Define input of the world tiff file
    filename_world = os.path.join(Dir_trash, "product", "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif")
    
    try:
        # Clip data to user extend
        data, Geo_out = RC.clip_data(filename_world, latlim, lonlim)
        
        # Save data of clipped array
        DC.Save_as_tiff(output_filename, data, Geo_out, 4326)
        
    except:
        
        RC.Clip_Dataset_GDAL(RC.clip_data(filename_world, latlim, lonlim))
    
    # Remove trash folder
    shutil.rmtree(Dir_trash)

    return()