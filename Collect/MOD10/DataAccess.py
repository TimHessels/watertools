# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MOD10
"""

# import general python modules
import os
import numpy as np
import pandas as pd
from osgeo import gdal
import urllib
from bs4 import BeautifulSoup
import re
import math
import datetime
import requests
import glob
from joblib import Parallel, delayed
import sys
if sys.version_info[0] == 3:
    import urllib.request    
    import urllib.parse
if sys.version_info[0] == 2:
    import urlparse

# Water Accounting modules
import watertools
import watertools.General.raster_conversions as RC
import watertools.General.data_conversions as DC

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar, cores, hdf_library, remove_hdf, period = "8-daily"):
    """
    This function downloads MOD10 8-daily data

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax] (values must be between -90 and 90)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    cores -- The number of cores used to run the routine. It can be 'False'
             to avoid using parallel computing routines.
	 nameDownload -- The name of the subset that must be download can be Fpar_500m or Lai_500m
    Waitbar -- 1 (Default) will print a waitbar
    """

    # Check start and end date and otherwise set the date to max
    if not Startdate:
        Startdate = pd.Timestamp('2000-02-18')
    if not Enddate:
        Enddate = pd.Timestamp('Now')


    # Make an array of the days of which the FPAR is taken
    if period == "8-daily":
        Dates = Make_TimeStamps(Startdate,Enddate)
    if period == "daily":
        Dates = pd.date_range(Startdate,Enddate,freq="D")
        
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

    # Make directory for the MODIS FPAR data
    Dir = Dir.replace("/", os.sep)
    if period == "8-daily":
        output_folder = os.path.join(Dir, 'MOD10')
    if period == "daily":    
        output_folder = os.path.join(Dir, 'MOD10', 'Daily')    

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # This function converts the values in the text file into horizontal and vertical number of the tiles which must be downloaded to cover the extent defined by the user
    TilesVertical, TilesHorizontal = watertools.Collect.MOD15.DataAccess.Get_tiles_from_txt(output_folder, hdf_library, latlim, lonlim)

    # Pass variables to parallel function and run
    args = [output_folder, TilesVertical, TilesHorizontal,lonlim, latlim, hdf_library, period]
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
        #   os.remove(os.path.join(output_folder, f))

    return(results)

def RetrieveData(Date, args):
    """
    This function retrieves MOD15 FPAR data for a given date from the
    http://e4ftl01.cr.usgs.gov/ server.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """
    # Argument
    [output_folder, TilesVertical, TilesHorizontal,lonlim, latlim, hdf_library, period] = args
    
    if period == "8-daily":
        FPARfileName = os.path.join(output_folder, 'SnowFrac_MOD10_-_8-daily_'  + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')
    if period == "daily":
        FPARfileName = os.path.join(output_folder, 'SnowFrac_MOD10_-_daily_'  + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')
 
    if not os.path.exists(FPARfileName):

        # Collect the data from the MODIS webpage and returns the data and lat and long in meters of those tiles
        try:
            Collect_data(TilesHorizontal, TilesVertical, Date, output_folder, hdf_library, period)
        except:
            print("Was not able to download the file")
        try:
            # Define the output name of the collect data function
            name_collect = os.path.join(output_folder, 'Merged.tif')
        
            # Reproject the MODIS product to epsg_to
            epsg_to ='4326'
            resolution = 0.005                
                    
            name_reprojected = RC.reproject_MODIS(name_collect, epsg_to, resolution = resolution)            

        
            # Clip the data to the users extend
            data, geo, proj = RC.clip_data(name_reprojected, latlim, lonlim)
        
            # Save the file as tiff
            DC.Save_as_tiff(name=FPARfileName, data=data, geo=geo, projection='WGS84')
        
            # remove the side products
            os.remove(os.path.join(output_folder, name_collect))
            os.remove(os.path.join(output_folder, name_reprojected))
                
        except:
            print("Failed for date: %s" %Date)
        
    return True

def Make_TimeStamps(Startdate,Enddate):
    '''
    This function determines all time steps of which the FPAR must be downloaded
    The time stamps are 8 daily.

    Keywords arguments:
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    '''

    # Define the DOY and year of the start day
    DOY = datetime.datetime.strptime(Startdate,'%Y-%m-%d').timetuple().tm_yday
    Year = datetime.datetime.strptime(Startdate,'%Y-%m-%d').timetuple().tm_year

    # Define the year of the end day
    YearEnd = datetime.datetime.strptime(Enddate,'%Y-%m-%d').timetuple().tm_year

    # Change the DOY of the start day into a DOY of MODIS day (8-daily) and create new startdate
    DOYstart = int(math.floor(DOY / 8.0) * 8) + 1
    DOYstart = str('%s-%s' %(DOYstart, Year))
    Day = datetime.datetime.strptime(DOYstart, '%j-%Y')
    Month = '%02d' % Day.month
    Day = '%02d' % Day.day
    Startdate = (str(Year) + '-' + str(Month) + '-' + str(Day))

    # Create the start and end data for the whole year
    YearStartDate = pd.date_range(Startdate, Enddate, freq = 'AS')
    YearEndDate = pd.date_range(Startdate, Enddate, freq = 'A')

    # Define the amount of years that are involved
    AmountOfYear = YearEnd - Year

    # If the startday is not in the same year as the enddate
    if AmountOfYear > 0:
        for i in range(0, AmountOfYear+1):
            if i == 0:
                Startdate1 = Startdate
                Enddate1 = YearEndDate[0]
                Dates = pd.date_range(Startdate1, Enddate1, freq = '8D')
            if i == AmountOfYear:
                Startdate1 = YearStartDate[-1]
                Enddate1 = Enddate
                Dates1 = pd.date_range(Startdate1, Enddate1, freq = '8D')
                Dates = Dates.union(Dates1)
            if (i != AmountOfYear and i != 0):
                Startdate1 = YearStartDate[i-AmountOfYear-1]
                Enddate1 = YearEndDate[i]
                Dates1 = pd.date_range(Startdate1, Enddate1, freq = '8D')
                Dates = Dates.union(Dates1)

    # If the startday is in the same year as the enddate
    if AmountOfYear == 0:
        Dates = pd.date_range(Startdate, Enddate, freq = '8D')

    return(Dates)


def Collect_data(TilesHorizontal,TilesVertical,Date,output_folder, hdf_library, period):
    '''
    This function downloads all the needed MODIS tiles from https://n5eil01u.ecs.nsidc.org/MOST/MOD10A2.006/ as a hdf file.

    Keywords arguments:
    TilesHorizontal -- [TileMin,TileMax] max and min horizontal tile number
    TilesVertical -- [TileMin,TileMax] max and min vertical tile number
    Date -- 'yyyy-mm-dd'
    output_folder -- 'C:/file/to/path/'
    '''

    # Make a new tile for the data
    sizeX = int((TilesHorizontal[1] - TilesHorizontal[0] + 1) * 2400)
    sizeY = int((TilesVertical[1] - TilesVertical[0] + 1) * 2400)
    DataTot = np.zeros((sizeY, sizeX))

    # Load accounts
    username, password = watertools.Functions.Random.Get_Username_PWD.GET('NASA')

    # Download the MODIS FPAR data
    if period == "8-daily":
        url = 'https://n5eil01u.ecs.nsidc.org/MOST/MOD10A2.061/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'
    if period == "daily":
        url = 'https://n5eil01u.ecs.nsidc.org/MOST/MOD10A1.061/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'
        
    dataset = requests.get(url, allow_redirects=False,stream = True)
    try:
        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True).content
    except:
        from requests.packages.urllib3.exceptions import InsecureRequestWarning
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
        get_dataset  = requests.get(dataset.headers['location'], auth = (username, password), verify = False).content

    soup = BeautifulSoup(get_dataset, "html.parser")

    if len(str(soup)) < 300:
        print('Download was not succesfull, please check NASA account')
        sys.exit(1)

    # Create the Lat and Long of the MODIS tile in meters
    for Vertical in range(int(TilesVertical[0]), int(TilesVertical[1])+1):
        Distance = 231.65635826395834*2 # resolution of a MODIS pixel in meter
        countY=(TilesVertical[1] - TilesVertical[0] + 1) - (Vertical - TilesVertical[0])

        for Horizontal in range(int(TilesHorizontal[0]), int(TilesHorizontal[1]) + 1):
            countX=Horizontal - TilesHorizontal[0] + 1

            for i in soup.findAll('a', attrs = {'href': re.compile('(?i)(hdf)$')}):

                # Find the file with the wanted tile number
                Vfile=str(i)[30:32]
                Hfile=str(i)[27:29]
                if int(Vfile) is int(Vertical) and int(Hfile) is int(Horizontal):

                    # Define the whole url name
                    if sys.version_info[0] == 3:
                        full_url = urllib.parse.urljoin(url, i['href'])

                    if sys.version_info[0] == 2:
                        full_url = urlparse.urljoin(url, i['href'])

		              # Reset the begin parameters for downloading
                    downloaded = 0
                    N=0

                    # if not downloaded try to download file
                    while downloaded == 0:

                        try:# open http and download whole .hdf
                            nameDownload_url = full_url
                            file_name = os.path.join(output_folder,nameDownload_url.split('/')[-1])
                            if os.path.isfile(file_name):
                                downloaded = 1
                            else:
                                x = requests.get(nameDownload_url, allow_redirects = False)
                                try:
                                    y = requests.get(x.headers['location'], auth = (username, password))
                                except:
                                    from requests.packages.urllib3.exceptions import InsecureRequestWarning
                                    requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

                                    y = requests.get(x.headers['location'], auth = (username, password), verify = False)
                                z = open(file_name, 'wb')
                                z.write(y.content)
                                z.close()
                                statinfo = os.stat(file_name)
                                # Say that download was succesfull
                                if int(statinfo.st_size) > 1000:
                                     downloaded = 1

                        # If download was not succesfull
                        except:

                            # Try another time
                            N = N + 1

						      # Stop trying after 10 times
                        if N == 10:
                            print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                            downloaded = 1

                    try:
                        # Open .hdf only band with SnowFrac and collect all tiles to one array
                        scale_factor = 1
                        dataset = gdal.Open(file_name)
                        sdsdict = dataset.GetMetadata('SUBDATASETS')
                        sdslist = [sdsdict[k] for k in sdsdict.keys() if '_1_NAME' in k]
                        sds = []

                        for n in sdslist:
                            sds.append(gdal.Open(n))
                            full_layer = [i for i in sdslist if 'MOD_Grid_Snow_500m' in i]

                            idx = sdslist.index(full_layer[0])
                            if Horizontal == TilesHorizontal[0] and Vertical == TilesVertical[0]:
                                geo_t = sds[idx].GetGeoTransform()

                                # get the projection value
                                proj = sds[idx].GetProjection()

                            data = sds[idx].ReadAsArray()
                            countYdata = (TilesVertical[1] - TilesVertical[0] + 2) - countY
                            DataTot[int((countYdata - 1) * 2400):int(countYdata * 2400), int((countX - 1) * 2400):int(countX * 2400)]=data * scale_factor
                        del data

                    # if the tile not exists or cannot be opened, create a nan array with the right projection
                    except:
                        if Horizontal==TilesHorizontal[0] and Vertical==TilesVertical[0]:
                             x1 = (TilesHorizontal[0] - 19) * 2400 * Distance
                             x4 = (TilesVertical[0] - 9) * 2400 * -1 * Distance
                             geo = 	[x1, Distance, 0.0, x4, 0.0, -Distance]
                             geo_t=tuple(geo)

                        proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
                        data=np.ones((2400, 2400)) * (-9999)
                        countYdata=(TilesVertical[1] - TilesVertical[0] + 2) - countY
                        DataTot[(countYdata - 1) * 2400:countYdata * 2400,(countX - 1) * 2400:countX * 2400] = data * 0.01


    # Make geotiff file
    name2 = os.path.join(output_folder, 'Merged.tif')
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name2, DataTot.shape[1], DataTot.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
    try:
        dst_ds.SetProjection(proj)
    except:
        proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
        x1 = (TilesHorizontal[0] - 18) * 2400 * Distance
        x4 = (TilesVertical[0] - 9) * 2400 * -1 * Distance
        geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
        geo_t = tuple(geo)
        dst_ds.SetProjection(proj)

    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo_t)
    dst_ds.GetRasterBand(1).WriteArray(DataTot)
    dst_ds = None
    sds = None
    return()

