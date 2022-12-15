# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MYD9
"""

# import general python modules
import os
import numpy as np
import pandas as pd
from osgeo import gdal
import urllib
from bs4 import BeautifulSoup
import re
import glob
import requests
from joblib import Parallel, delayed
import sys
if sys.version_info[0] == 3:
    import urllib.request    
    import urllib.parse
if sys.version_info[0] == 2:
    import urlparse
    import urllib2


# Water Accounting modules
import watertools
import watertools.General.raster_conversions as RC
import watertools.General.data_conversions as DC

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar, band, resolution, cores, hdf_library, remove_hdf):
    """
    This function downloads MOD9 daily data

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax] (values must be between -90 and 90)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    cores -- The number of cores used to run the routine. It can be 'False'
             to avoid using parallel computing routines.
    Waitbar -- 1 (Default) will print a waitbar
    """

    # Check start and end date and otherwise set the date to max
    if not Startdate:
        Startdate = pd.Timestamp('2000-02-24')
    if not Enddate:
        Enddate = pd.Timestamp('Now')

    # Make an array of the days of which the NDVI is taken
    Dates = pd.date_range(Startdate, Enddate, freq = 'D')

    # Create Waitbar
    if Waitbar == 1:
        import watertools.Functions.Random.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -90 or latlim[1] > 90:
        print('Latitude above 90N or below 90S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -90)
        latlim[1] = np.min(latlim[1], 90)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)

    # Make directory for the MODIS NDVI data
    Dir = Dir.replace("/", os.sep)
    output_folder = os.path.join(Dir, 'Reflectance', 'MYD9', 'Band_%s_%s' %(band,resolution))
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    TilesVertical, TilesHorizontal = watertools.Collect.MOD15.DataAccess.Get_tiles_from_txt(output_folder, hdf_library, latlim, lonlim)

    # Pass variables to parallel function and run
    args = [output_folder, TilesVertical, TilesHorizontal, lonlim, latlim, band, resolution, hdf_library]
    if not cores:
        for Date in Dates:
            RetrieveData(Date, args)
            if Waitbar == 1:
                amount += 1
                WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
        results = True
    else:
        results = Parallel(n_jobs=cores)(delayed(RetrieveData)(Date, args)
                                         for Date in Dates)
    if remove_hdf == 1:
        # Remove all .hdf files
        os.chdir(output_folder)
        files = glob.glob("*.hdf")
        for f in files:
            os.remove(os.path.join(output_folder, f))

        # Remove all .txt files
        #files = glob.glob("*.txt")
        #for f in files:
        #    os.remove(os.path.join(output_folder, f))

    return(results)

def RetrieveData(Date, args):
    """
    This function retrieves MOD9 Reflectance data for a given date from the
    http://e4ftl01.cr.usgs.gov/ server.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """
    # Argument
    [output_folder, TilesVertical, TilesHorizontal, lonlim, latlim, band, resolution, hdf_library] = args

    if resolution == "250m":
        letter = "Q"
    else:                
        letter = "A"
        
    if band == "state":
        ReffileName = os.path.join(output_folder, 'ReflectanceBand_%s_MYD09G%s_-_daily_'%(band, letter) + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')        
    else:
        ReffileName = os.path.join(output_folder, 'ReflectanceBand%d_MYD09G%s_-_daily_'%(band, letter) + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')

    if not os.path.exists(ReffileName):    
        # Collect the data from the MODIS webpage and returns the data and lat and long in meters of those tiles
        try:
            Collect_data(TilesHorizontal, TilesVertical, Date, output_folder, band, resolution, hdf_library)
        except:
            print("Was not able to download the file")
    
        # Define the output name of the collect data function
        name_collect = os.path.join(output_folder, 'Merged.tif')
        try:
            # Reproject the MODIS product to epsg_to
            epsg_to ='4326'

            if resolution == "250m":
                resolution = 0.0025
            elif resolution == "1000m":
                resolution = 0.01 
            else:
                resolution = 0.005               
                    
            name_reprojected = RC.reproject_MODIS(name_collect, epsg_to, resolution = resolution)
        
            # Clip the data to the users extend
            data, geo, proj = RC.clip_data(name_reprojected, latlim, lonlim)
        
            # Save results as Gtiff
            DC.Save_as_tiff(name=ReffileName, data=data, geo=geo, projection='WGS84')
        
            # remove the side products
            os.remove(os.path.join(output_folder, name_collect))
            os.remove(os.path.join(output_folder, name_reprojected))
        except:
           print('data for %02d-%02d-%d is not available' %(Date.day, Date.month, Date.year))
            
    return True


def Collect_data(TilesHorizontal,TilesVertical,Date,output_folder, band, resolution, hdf_library):
    '''
    This function downloads all the needed MODIS tiles from http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.006/ as a hdf file.

    Keywords arguments:
    TilesHorizontal -- [TileMin,TileMax] max and min horizontal tile number
    TilesVertical -- [TileMin,TileMax] max and min vertical tile number
    Date -- 'yyyy-mm-dd'
    output_folder -- 'C:/file/to/path/'
    '''
    
    if band =="state":
        if resolution != "1000m":
            sys.exit('Bands State, are only available in 1000m resolution')
            
        if resolution == "1000m":
            size_factor = 4

    else:
        if band>=3 and resolution == "250m":
            sys.exit('Bands higher than 3, are only available in 500m resolution')
            
        if resolution == "250m":
            size_factor = 1
        else:
            size_factor = 2    
        
    # Make a new tile for the data
    sizeX = int((TilesHorizontal[1] - TilesHorizontal[0] + 1) * 4800/size_factor)
    sizeY = int((TilesVertical[1] - TilesVertical[0] + 1) * 4800/size_factor)
    DataTot = np.zeros((sizeY, sizeX))

    # Load accounts
    BEARER = watertools.Functions.Random.Get_Username_PWD.GET('NASA_BEARER')

    s = requests.Session()  
    headers = {
        'Authorization': 'Bearer %s'%BEARER[0]
        }                

    # Download the MODIS NDVI data
    if resolution == "250m":
        #url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD09GQ.061/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'
        url = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MYD09GQ/%s/%03s" %(Date.strftime('%Y'),  Date.strftime('%j'))
        letter = "Q"
    else:                
        #url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD09GA.061/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'
        url = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MYD09GA/%s/%03s" %(Date.strftime('%Y'),  Date.strftime('%j'))
        letter = "A"

    # Create the Lat and Long of the MODIS tile in meters
    for Vertical in range(int(TilesVertical[0]), int(TilesVertical[1])+1):
        Distance = 231.65635826395834 * size_factor # resolution of a MODIS pixel in meter
        countY=(TilesVertical[1] - TilesVertical[0] + 1) - (Vertical - TilesVertical[0])

        for Horizontal in range(int(TilesHorizontal[0]), int(TilesHorizontal[1]) + 1):
            countX=Horizontal - TilesHorizontal[0] + 1

		    # Reset the begin parameters for downloading
            downloaded = 0
            N=0

	         # Check the library given by user
            if hdf_library is not None:
                os.chdir(hdf_library)
                hdf_name = glob.glob("MYD09G%s.A%s%03s.h%02dv%02d.*" %(letter, Date.strftime('%Y'), Date.strftime('%j'), Horizontal, Vertical))

                if len(hdf_name) == 1:
                    hdf_file = os.path.join(hdf_library, hdf_name[0])

                    if os.path.exists(hdf_file):
                        downloaded = 1
                        file_name = hdf_file

            if not downloaded == 1:
                
                r = s.get(''.join([url,'?fields=all&format=json']), headers=headers, timeout = 2000)
                if r.status_code == 200:            
                    f = r.text 
                else:
                    print("Got an ERROR 404 not connected!!!")
                    
                # Sum all the files on the server
                soup = BeautifulSoup(f, "lxml")    
                    
                try:

                    for i in soup.findAll('a', attrs = {'href': re.compile('(?i)(hdf)$')}):
    
                        # Find the file with the wanted tile number
                        Vfile=str(i)[80:82]
                        Hfile=str(i)[77:79]
                        if int(Vfile) is int(Vertical) and int(Hfile) is int(Horizontal):
                       
                            # Define the whole url name
                            if sys.version_info[0] == 3:
                                full_url = urllib.parse.urljoin(url, i['href'])
    
                            if sys.version_info[0] == 2:
                                full_url = urlparse.urljoin(url, i['href'])
    
                            # if not downloaded try to download file
                            while downloaded == 0:
    
                                try:# open http and download whole .hdf
                                    nameDownload = full_url
                                    file_name = os.path.join(output_folder, nameDownload.split('/')[-1])
                                    if os.path.isfile(file_name):
                                        downloaded = 1
                                    else:
                                                           
                                        r = s.get(nameDownload, headers=headers, stream=True, timeout = 2000)
                                        with open(file_name, 'wb') as f:
                                            for chunk in r.iter_content(chunk_size=1024 * 1024):
                                                if chunk:  # filter out keep-alive new chunks
                                                    f.write(chunk)
                                                    # f.flush()                                       
                                                              
                                        statinfo = os.stat(file_name)
                                        # Say that download was succesfull
                                        if int(statinfo.st_size) > 10000:
                                             downloaded = 1
                                        else:
                                            print("Size is too small")
    
                                # If download was not succesfull
                                except:
    
                                    # Try another time
                                    N = N + 1
    
        				        # Stop trying after 10 times
                                if N == 10:
                                    print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                                    downloaded = 1

                except:
                        print("Url not found: %s" %url)                                                 
                                
            try:
                # Open .hdf only band with NDVI and collect all tiles to one array
                dataset = gdal.Open(file_name)
                sdsdict = dataset.GetMetadata('SUBDATASETS')
                if resolution =="250m":
                    sdslist = [sdsdict[k] for k in sdsdict.keys() if 'SUBDATASET_%d_NAME'%int(band+1) in k]  
                elif resolution == "1000m":
                    sdslist = [sdsdict[k] for k in sdsdict.keys() if 'SUBDATASET_2_NAME' in k]  
                else:
                    sdslist = [sdsdict[k] for k in sdsdict.keys() if 'SUBDATASET_%d_NAME'%int(band+11) in k]
                sds = []

                for n in sdslist:
                    sds.append(gdal.Open(n))
                    if resolution != "1000m":
                        full_layer = [i for i in sdslist if 'sur_refl_b%02d_1'%band in i]     
                    else:
                        full_layer = [i for i in sdslist if 'state_1km_1' in i]     
                        
                    idx = sdslist.index(full_layer[0])
                    if Horizontal == TilesHorizontal[0] and Vertical == TilesVertical[0]:
                        geo_t = sds[idx].GetGeoTransform()

                        # get the projection value
                        proj = sds[idx].GetProjection()

                    data = sds[idx].ReadAsArray()
                    countYdata = (TilesVertical[1] - TilesVertical[0] + 2) - countY
                    DataTot[int((countYdata - 1) * 4800/size_factor):int(countYdata * 4800/size_factor), int((countX - 1) * 4800/size_factor):int(countX * 4800/size_factor)]=data
                del data

            # if the tile not exists or cannot be opened, create a nan array with the right projection
            except:
                if Horizontal==TilesHorizontal[0] and Vertical==TilesVertical[0]:
                     x1 = (TilesHorizontal[0] - 19) * 4800/size_factor * Distance
                     x4 = (TilesVertical[0] - 9) * 4800/size_factor* -1 * Distance
                     geo = 	[x1, Distance, 0.0, x4, 0.0, -Distance]
                     geo_t=tuple(geo)

                proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
                data=np.ones((int(4800/size_factor),int(4800/size_factor))) * (-9999)
                countYdata=(TilesVertical[1] - TilesVertical[0] + 2) - countY
                DataTot[int((countYdata - 1) * 4800/size_factor):int(countYdata * 4800/size_factor),int((countX - 1) * 4800/size_factor):int(countX * 4800/size_factor)] = data

    # Make geotiff file
    name2 = os.path.join(output_folder, 'Merged.tif')
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name2, DataTot.shape[1], DataTot.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
    try:
         dst_ds.SetProjection(proj)
    except:
        proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
        x1 = (TilesHorizontal[0] - 18) * 4800/size_factor * Distance
        x4 = (TilesVertical[0] - 9) * 4800/size_factor * -1 * Distance
        geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
        geo_t = tuple(geo)
        dst_ds.SetProjection(proj)

    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo_t)
    dst_ds.GetRasterBand(1).WriteArray(DataTot*0.0001)
    dst_ds = None
    sds = None
    return()
