# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MOD11
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
import glob
import datetime
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
import watertools.General.raster_conversions as RC
import watertools.General.data_conversions as DC
import watertools

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, TimeStep, Waitbar, cores, hdf_library, remove_hdf, angle_info = 0, time_info = 0, day_night = "day"):
    """
    This function downloads MOD13 16-daily data

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
        Startdate = pd.Timestamp('2000-02-18')
    if not Enddate:
        Enddate = pd.Timestamp('Now')

    # Make an array of the days of which the LST is taken
    if TimeStep == 8:
        Dates = Make_TimeStamps(Startdate,Enddate)
        if day_night == "night":
            TimeStepName = "8_Night"        
        else:
            TimeStepName = '8_Daily'
    if TimeStep == 1:
        Dates = pd.date_range(Startdate,Enddate,freq = 'D')
        if day_night == "night":
            TimeStepName = "Night"
        else:
            TimeStepName = 'Daily'

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
    output_folder = os.path.join(Dir, 'LST', 'MOD11', '%s' %TimeStepName)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define which MODIS tiles are required
    TilesVertical, TilesHorizontal = watertools.Collect.MOD15.DataAccess.Get_tiles_from_txt(output_folder, hdf_library, latlim, lonlim)

    # Pass variables to parallel function and run
    args = [output_folder, TilesVertical, TilesHorizontal,lonlim, latlim, TimeStep, hdf_library, angle_info, time_info, day_night]
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
    This function retrieves MOD11 LST data for a given date from the
    https://e4ftl01.cr.usgs.gov/ server.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """
    # Argument
    [output_folder, TilesVertical, TilesHorizontal,lonlim, latlim, TimeStep, hdf_library, angle_info, time_info, day_night] = args

    if TimeStep == 8:
        LSTfileNamePart = os.path.join(output_folder, 'LST_MOD11A2_K_8-daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')
    if TimeStep == 1:
        LSTfileNamePart = os.path.join(output_folder, 'LST_MOD11A1_K_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.*.tif')
    if angle_info == 1:
        OnsangfileNamePart = os.path.join(output_folder, 'Angle_MOD11A1_degrees_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')    
    if time_info == 1:
        TimefileNamePart = os.path.join(output_folder, 'Time_MOD11A1_hour_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')    
  
    filesMOD = glob.glob(LSTfileNamePart)
    if angle_info == 1:
        filesANGLE = glob.glob(OnsangfileNamePart)
    else:
        filesANGLE = ["not_required"]
        
    if time_info == 1:
        filesTime = glob.glob(TimefileNamePart)
    else:
        filesTime = ["not_required"]    
        
    if not (len(filesMOD) == 1 and len(filesANGLE) == 1 and len(filesTime) == 1):

        # Collect the data from the MODIS webpage and returns the data and lat and long in meters of those tiles
        try:
            Collect_data(TilesHorizontal, TilesVertical, Date, output_folder, TimeStep, hdf_library, day_night, angle_info)
        except:
            print("Was not able to download the file")
        try:
            # Define the output name of the collect data function
            name_collect = os.path.join(output_folder, 'Merged.tif')
        
            # Reproject the MODIS product to epsg_to
            epsg_to ='4326'
            resolution = 0.01
            name_reprojected = RC.reproject_MODIS(name_collect, epsg_to, resolution = resolution)
        
            # Clip the data to the users extend
            data, geo, proj = RC.clip_data(name_reprojected, latlim, lonlim)
                
            # Save results as Gtiff
            if TimeStep == 8:
                LSTfileName = os.path.join(output_folder, 'LST_MOD11A2_K_8-daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')
            if TimeStep == 1:
                name_collect_time = os.path.join(output_folder, 'Merged_Time.tif')
                name_reprojected_time = RC.reproject_MODIS(name_collect_time, epsg_to, resolution = resolution)
                data_time, geo, proj = RC.clip_data(name_reprojected_time, latlim, lonlim)
                data_time[data_time==25.5] = np.nan
                data_time_ave = np.nanmean(data_time)
                try:
                    hour_GMT = int(np.floor(data_time_ave))
                    minutes_GMT = int((data_time_ave - np.floor(data_time_ave))*60)    
                except:
                    hour_GMT = int(10)
                    minutes_GMT = int(30)
                LSTfileName = os.path.join(output_folder, 'LST_MOD11A1_K_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.%02d%02d.tif'%(hour_GMT,minutes_GMT))
                os.remove(name_collect_time)
                os.remove(name_reprojected_time) 
                if time_info == 1: 
                    TimefileName = os.path.join(output_folder, 'Time_MOD11A1_hour_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')    
                    DC.Save_as_tiff(name=TimefileName, data=data_time, geo=geo, projection='WGS84')

                if angle_info == 1:
                    name_collect_angle = os.path.join(output_folder, 'Merged_Obsang.tif')
                    name_reprojected_angle = RC.reproject_MODIS(name_collect_angle, epsg_to, resolution = resolution)
                    data_angle, geo, proj = RC.clip_data(name_reprojected_angle, latlim, lonlim)
                    data_angle[data_angle==25.5] = np.nan
                    OnsangfileName = os.path.join(output_folder, 'Angle_MOD11A1_degrees_daily_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')    
                    data_angle[data_angle==0.] = -9999
                    DC.Save_as_tiff(name=OnsangfileName, data=data_angle, geo=geo, projection='WGS84')
                    os.remove(name_collect_angle)
                    os.remove(name_reprojected_angle) 
                    
            data[data==0.] = -9999
            DC.Save_as_tiff(name=LSTfileName, data=data, geo=geo, projection='WGS84')
        
            # remove the side products
            os.remove(os.path.join(output_folder, name_collect))
            os.remove(os.path.join(output_folder, name_reprojected))
            
        except Exception as e:
            print("Failed for date: %s" %Date)
            print(e)
            
    return True

def Make_TimeStamps(Startdate,Enddate):
    '''
    This function determines all time steps of which the LST must be downloaded
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

    # Change the DOY of the start day into a DOY of MODIS day (16-daily) and create new startdate
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

def Collect_data(TilesHorizontal,TilesVertical,Date,output_folder, TimeStep, hdf_library, day_night, angle_info = 0):
    '''
    This function downloads all the needed MODIS tiles from http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.006/ as a hdf file.

    Keywords arguments:
    TilesHorizontal -- [TileMin,TileMax] max and min horizontal tile number
    TilesVertical -- [TileMin,TileMax] max and min vertical tile number
    Date -- 'yyyy-mm-dd'
    output_folder -- 'C:/file/to/path/'
    '''

    # Make a new tile for the data
    sizeX = int((TilesHorizontal[1] - TilesHorizontal[0] + 1) * 1200)
    sizeY = int((TilesVertical[1] - TilesVertical[0] + 1) * 1200)
    DataTot = np.zeros((sizeY, sizeX))
    if angle_info == 1:
        DataTot_ObsAng = np.zeros((sizeY, sizeX))        
    if TimeStep == 1:    
        DataTot_Time = np.zeros((sizeY, sizeX))
        
    # Load accounts
    username, password = watertools.Functions.Random.Get_Username_PWD.GET('NASA')

    # Create the Lat and Long of the MODIS tile in meters
    for Vertical in range(int(TilesVertical[0]), int(TilesVertical[1])+1):
        Distance = 4*231.65635826395834 # resolution of a MODIS pixel in meter
        countY=int((TilesVertical[1] - TilesVertical[0] + 1) - (Vertical - TilesVertical[0]))

        for Horizontal in range(int(TilesHorizontal[0]), int(TilesHorizontal[1]) + 1):
            countX=int(Horizontal - TilesHorizontal[0] + 1)

            # Download the MODIS NDVI data
            if TimeStep == 8:
                url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.061/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'
            if TimeStep == 1:
                url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.061/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'

		    # Reset the begin parameters for downloading
            downloaded = 0
            N=0

			   # Check the library given by user
            if hdf_library is not None:
                os.chdir(hdf_library)
                if TimeStep == 8:
                    hdf_name = glob.glob("MOD11A2.A%s%03s.h%02dv%02d.*" %(Date.strftime('%Y'), Date.strftime('%j'), Horizontal, Vertical))
                if TimeStep == 1:
                    hdf_name = glob.glob("MOD11A1.A%s%03s.h%02dv%02d.*" %(Date.strftime('%Y'), Date.strftime('%j'), Horizontal, Vertical))

                if len(hdf_name) == 1:
                    hdf_file = os.path.join(hdf_library, hdf_name[0])

                    if os.path.exists(hdf_file):
                        downloaded = 1
                        file_name = hdf_file

            if not downloaded == 1:

                try:
                    # Get files on FTP server
                    if sys.version_info[0] == 3:
                        f = urllib.request.urlopen(url)
    
                    if sys.version_info[0] == 2:
                        f = urllib2.urlopen(url)

                    # Sum all the files on the server
                    soup = BeautifulSoup(f, "lxml")
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
    
                            # if not downloaded try to download file
                            while downloaded == 0:
    
                                try:# open http and download whole .hdf
                                    nameDownload = full_url
                                    file_name = os.path.join(output_folder,nameDownload.split('/')[-1])
                                    if os.path.isfile(file_name):
                                        print("file ", file_name, " already exists")
                                        downloaded = 1
                                    else:
                                        x = requests.get(nameDownload, allow_redirects = False)
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
                                        if int(statinfo.st_size) > 10000:
                                             downloaded = 1
    
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
                if day_night == "day":
                    sdslist = [sdsdict[k] for k in sdsdict.keys() if (('SUBDATASET_1_NAME') in k or ('SUBDATASET_3_NAME') in k or ('SUBDATASET_4_NAME') in k)]
                else:
                    sdslist = [sdsdict[k] for k in sdsdict.keys() if (('SUBDATASET_5_NAME') in k or ('SUBDATASET_7_NAME') in k or ('SUBDATASET_8_NAME') in k)]                    
                sds = []
                sds_time = []
                sds_obsang = []
                
                for n in sdslist:
                    sds.append(gdal.Open(n))
                    if day_night == "day":
                        full_layer = [i for i in sdslist if 'LST_Day_1km' in i]
                    else:
                        full_layer = [i for i in sdslist if 'LST_Night_1km' in i]                       
                    idx = sdslist.index(full_layer[0])
                    if Horizontal == TilesHorizontal[0] and Vertical == TilesVertical[0]:
                        geo_t = sds[idx].GetGeoTransform()

                        # get the projection value
                        proj = sds[idx].GetProjection()

                    data = sds[idx].ReadAsArray()
                    countYdata = int((TilesVertical[1] - TilesVertical[0] + 2) - countY)
                    DataTot[int((countYdata - 1) * 1200):int(countYdata * 1200), int((countX - 1) * 1200):int(countX * 1200)]=data * 0.02
                del data
 
                if TimeStep == 1:
                    if day_night == "day":
                        full_layer_time = [i for i in sdslist if 'Day_view_time' in i]
                    else:
                        full_layer_time = [i for i in sdslist if 'Night_view_time' in i]                        
                    idx_time = sdslist.index(full_layer_time[0])
                    sds_time.append(gdal.Open(sdslist[idx_time]))                
                    data_time = sds_time[0].ReadAsArray()
                    DataTot_Time[int((countYdata - 1) * 1200):int(countYdata * 1200), int((countX - 1) * 1200):int(countX * 1200)]=data_time * 0.1
                    del data_time

                if angle_info == 1:
                    if day_night == "day":
                        full_layer_obsang = [i for i in sdslist if 'Day_view_angl' in i]
                    else:
                        full_layer_obsang = [i for i in sdslist if 'Night_view_angl' in i]                        
                    idx_obsang = sdslist.index(full_layer_obsang[0])
                    sds_obsang.append(gdal.Open(sdslist[idx_obsang]))            
                    data_obsang = sds_obsang[0].ReadAsArray()
                    DataTot_ObsAng[int((countYdata - 1) * 1200):int(countYdata * 1200), int((countX - 1) * 1200):int(countX * 1200)]= data_obsang                 
                    del data_obsang


            # if the tile not exists or cannot be opened, create a nan array with the right projection
            except:
                if Horizontal==TilesHorizontal[0] and Vertical==TilesVertical[0]:
                     x1 = (TilesHorizontal[0] - 19) * 1200 * Distance
                     x4 = (TilesVertical[0] - 9) * 1200 * -1 * Distance
                     geo = 	[x1, Distance, 0.0, x4, 0.0, -Distance]
                     geo_t=tuple(geo)

                proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
                data=np.ones((1200,1200)) * (-9999/0.02)
                countYdata=(TilesVertical[1] - TilesVertical[0] + 2) - countY
                DataTot[int((countYdata - 1) * 1200):int(countYdata * 1200),int((countX - 1) * 1200):int(countX * 1200)] = data * 0.02
                DataTot[DataTot < 1] = -9999
                
                if TimeStep == 1:
                    data_time=np.ones((1200,1200)) * (-9999/0.1)
                    DataTot_Time[int((countYdata - 1) * 1200):int(countYdata * 1200),int((countX - 1) * 1200):int(countX * 1200)] = data_time * 0.1

                if angle_info == 1:
                    data_obsang=np.ones((1200,1200)) * (-9999)
                    DataTot_ObsAng[int((countYdata - 1) * 1200):int(countYdata * 1200),int((countX - 1) * 1200):int(countX * 1200)] = data_obsang

    # set limits
    if angle_info == 1:
        DataTot_ObsAng[DataTot_ObsAng<255] = DataTot_ObsAng[DataTot_ObsAng<255] - 65
        DataTot_ObsAng[DataTot_ObsAng==255] = -9999
        DataTot_ObsAng[DataTot_ObsAng>255] = DataTot_ObsAng[DataTot_ObsAng>255] - 256


    # Make geotiff file
    DataTot[DataTot < 1] = -9999    
    name2 = os.path.join(output_folder, 'Merged.tif')
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name2, DataTot.shape[1], DataTot.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
    try:
         dst_ds.SetProjection(proj)
    except:
        proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
        x1 = (TilesHorizontal[0] - 18) * 1200 * Distance
        x4 = (TilesVertical[0] - 9) * 1200 * -1 * Distance
        geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
        geo_t = tuple(geo)
        dst_ds.SetProjection(proj)

    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo_t)
    dst_ds.GetRasterBand(1).WriteArray(DataTot)
    dst_ds = None
    sds = None

    if TimeStep == 1:    
        DataTot_Time[np.logical_and(DataTot_Time < 0., DataTot_Time>24.)] = -9999
        # Make geotiff file of the time data
        name_out_time = os.path.join(output_folder, 'Merged_Time.tif')
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(name_out_time, DataTot_Time.shape[1], DataTot_Time.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
        try:
             dst_ds.SetProjection(proj)
        except:
            proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
            x1 = (TilesHorizontal[0] - 18) * 1200 * Distance
            x4 = (TilesVertical[0] - 9) * 1200 * -1 * Distance
            geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
            geo_t = tuple(geo)
            dst_ds.SetProjection(proj)
    
        dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
        dst_ds.SetGeoTransform(geo_t)
        dst_ds.GetRasterBand(1).WriteArray(DataTot_Time)
        dst_ds = None
        sds = None        

    if angle_info == 1:
        name_out_obsang = os.path.join(output_folder, 'Merged_Obsang.tif')
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(name_out_obsang, DataTot_ObsAng.shape[1], DataTot_ObsAng.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
        try:
             dst_ds.SetProjection(proj)
        except:
            proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
            x1 = (TilesHorizontal[0] - 18) * 1200 * Distance
            x4 = (TilesVertical[0] - 9) * 1200 * -1 * Distance
            geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
            geo_t = tuple(geo)
            dst_ds.SetProjection(proj)
    
        dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
        dst_ds.SetGeoTransform(geo_t)
        dst_ds.GetRasterBand(1).WriteArray(DataTot_ObsAng)
        dst_ds = None
        sds = None    

    return()
