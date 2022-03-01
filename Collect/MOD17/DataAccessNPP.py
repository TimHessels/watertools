# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MOD17
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

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar, cores, hdf_library, remove_hdf):
    """
    This function downloads MOD17 yearly NPP data

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

    # Make an array of the days of which the NPP is taken
    yearstart = pd.Timestamp(Startdate).year
    yearend = pd.Timestamp(Enddate).year
    Startdate_NPP='%s-01-01' % yearstart
    Enddate_NPP='%s-12-31'% yearend
    Dates = pd.date_range(Startdate_NPP, Enddate_NPP, freq = 'AS')

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

    # Make directory for the MODIS NPP data
    Dir = Dir.replace("/", os.sep)
    output_folder = os.path.join(Dir, 'NPP', 'MOD17')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define which MODIS tiles are required
    TilesVertical, TilesHorizontal = watertools.Collect.MOD15.DataAccess.Get_tiles_from_txt(output_folder, hdf_library, latlim, lonlim)

    # Pass variables to parallel function and run
    args = [output_folder, TilesVertical, TilesHorizontal, lonlim, latlim, hdf_library]
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
    This function retrieves MOD17 NPP data for a given date from the
    http://e4ftl01.cr.usgs.gov/ server.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """
    # Argument
    [output_folder, TilesVertical, TilesHorizontal, lonlim, latlim, hdf_library] = args
    NPPfileName = os.path.join(output_folder, 'NPP_MOD17_kg-C-m^-2_yearly_' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '.tif')

    if not os.path.exists(NPPfileName):
            
        # Collect the data from the MODIS webpage and returns the data and lat and long in meters of those tiles
        try:
            Collect_data(TilesHorizontal, TilesVertical, Date, output_folder, hdf_library)
        except:
            print("Was not able to download the file")
        try:
            
            # Define the output name of the collect data function
            name_collect = os.path.join(output_folder, 'Merged.tif')
        
            # Reproject the MODIS product to epsg_to
            epsg_to ='4326'
            name_reprojected = RC.reproject_MODIS(name_collect, epsg_to, resolution = 0.005)
        
            # Clip the data to the users extend
            data, geo, proj = RC.clip_data(name_reprojected, latlim, lonlim)
        
            # Save results as Gtiff
            DC.Save_as_tiff(name=NPPfileName, data=data, geo=geo, projection='WGS84')
        
            # remove the side products
            os.remove(os.path.join(output_folder, name_collect))
            os.remove(os.path.join(output_folder, name_reprojected))
        except:
            print("Failed for date: %s" %Date)
            
    return True

def Collect_data(TilesHorizontal,TilesVertical,Date,output_folder, hdf_library):
    '''
    This function downloads all the needed MODIS tiles from http://e4ftl01.cr.usgs.gov/MOLT/MOD17A3H.006/ as a hdf file.

    Keywords arguments:
    TilesHorizontal -- [TileMin,TileMax] max and min horizontal tile number
    TilesVertical -- [TileMin,TileMax] max and min vertical tile number
    Date -- 'yyyy-mm-dd'
    output_folder -- 'C:/file/to/path/'
    '''

    # Make a new tile for the data
    # NPP_SIZE = 2
    NPP_SIZE = 4
    sizeX = int((TilesHorizontal[1] - TilesHorizontal[0] + 1) * 4800 / NPP_SIZE )
    sizeY = int((TilesVertical[1] - TilesVertical[0] + 1) * 4800 / NPP_SIZE)
    DataTot = np.zeros((sizeY, sizeX))

    # Load accounts
    username, password = watertools.Functions.Random.Get_Username_PWD.GET('NASA')

    # Create the Lat and Long of the MODIS tile in meters
    for Vertical in range(int(TilesVertical[0]), int(TilesVertical[1])+1):
        Distance = 231.65635826395834 * NPP_SIZE # resolution of a MODIS pixel in meter
        countY=(TilesVertical[1] - TilesVertical[0] + 1) - (Vertical - TilesVertical[0])

        for Horizontal in range(int(TilesHorizontal[0]), int(TilesHorizontal[1]) + 1):
            countX=Horizontal - TilesHorizontal[0] + 1

            # Download the MODIS NPP data
            #url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD17A3H.006/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'
            #url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD17A3.055/' + Date.strftime('%Y') + '.' + Date.strftime('%m') + '.' + Date.strftime('%d') + '/'
            url = 'http://files.ntsg.umt.edu/data/NTSG_Products/MOD17/MOD17A3/Y' + Date.strftime('%Y') + '/'

		      # Reset the begin parameters for downloading
            downloaded = 0
            N=0

	         # Check the library given by user
            if hdf_library is not None:
                os.chdir(hdf_library)
                hdf_name = glob.glob("MOD17A3.A%s%03s.h%02dv%02d.*" %(Date.strftime('%Y'), Date.strftime('%j'), Horizontal, Vertical))

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
                        Vfile=str(i)[-28:-26]
                        Hfile=str(i)[-31:-29]
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
                                        """
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
                                        """
                                        if sys.version_info[0] == 3:
                                            urllib.request.urlretrieve(nameDownload, file_name)
                                        if sys.version_info[0] == 2:
                                            urllib.urlretrieve(nameDownload, file_name)
    
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
                # Open .hdf only band with NPP and collect all tiles to one array
                dataset = gdal.Open(file_name)
                sdsdict = dataset.GetMetadata('SUBDATASETS')
                #sdslist = [sdsdict[k] for k in sdsdict.keys() if '_2_NAME' in k]
                sdslist = [sdsdict[k] for k in sdsdict.keys() if '_2_NAME' in k]
                sds = []

                for n in sdslist:
                    sds.append(gdal.Open(n))
                    #yfull_layer = [i for i in sdslist if 'Npp_500m' in i]
                    full_layer = [i for i in sdslist if 'Npp_1km' in i]
                    idx = sdslist.index(full_layer[0])
                    if Horizontal == TilesHorizontal[0] and Vertical == TilesVertical[0]:
                        geo_t = sds[idx].GetGeoTransform()

                        # get the projection value
                        proj = sds[idx].GetProjection()

                    data = sds[idx].ReadAsArray()
                    countYdata = (TilesVertical[1] - TilesVertical[0] + 2) - countY
                    DataTot[int((countYdata - 1) * 4800 / NPP_SIZE):int(countYdata * 4800 / NPP_SIZE), int((countX - 1) * 4800 / NPP_SIZE):int(countX * 4800 / NPP_SIZE)]=data * 0.0001
                del data

            # if the tile not exists or cannot be opened, create a nan array with the right projection
            except:
                if Horizontal==TilesHorizontal[0] and Vertical==TilesVertical[0]:
                     x1 = (TilesHorizontal[0] - 19) * 4800 / NPP_SIZE * Distance
                     x4 = (TilesVertical[0] - 9) * 4800 / NPP_SIZE * -1 * Distance
                     geo = 	[x1, Distance, 0.0, x4, 0.0, -Distance]
                     geo_t=tuple(geo)

                proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
                data=np.ones((int(4800 / NPP_SIZE), int(4800 / NPP_SIZE))) * (0)
                countYdata=(TilesVertical[1] - TilesVertical[0] + 2) - countY
                DataTot[int((countYdata - 1) * 4800 / NPP_SIZE):int(countYdata * 4800 / NPP_SIZE),int((countX - 1) * 4800 / NPP_SIZE):int(countX * 4800 / NPP_SIZE)] = data * 0.0001

    # Make geotiff file
    DataTot[DataTot>3.27]=-9999
    name2 = os.path.join(output_folder, 'Merged.tif')
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(name2, DataTot.shape[1], DataTot.shape[0], 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
    try:
         dst_ds.SetProjection(proj)
    except:
        proj='PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
        x1 = (TilesHorizontal[0] - 18) * 4800 / NPP_SIZE * Distance
        x4 = (TilesVertical[0] - 9) * 4800 / NPP_SIZE * -1 * Distance
        geo = [x1, Distance, 0.0, x4, 0.0, -Distance]
        geo_t = tuple(geo)
        dst_ds.SetProjection(proj)

    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.SetGeoTransform(geo_t)
    dst_ds.GetRasterBand(1).WriteArray(DataTot)
    dst_ds = None
    sds = None
    return()
