# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 19:04:22 2016

@author: tih
"""
import pandas as pd
import glob
import gdal
import osr
import os
import sys
import numpy as np
import subprocess
from pyproj import Proj, transform
import scipy.interpolate
import fiona
import requests
import pycurl
import shapefile

def Run_command_window(argument):
    """
    This function runs the argument in the command window without showing cmd window

    Keyword Arguments:
    argument -- string, name of the adf file
    """
    if os.name == 'posix':
        argument = argument.replace(".exe","")
        argument = argument.replace(".bat","")
        os.system(argument)

    else:
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW

        process = subprocess.Popen(argument, startupinfo=startupinfo, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        process.wait()

    return()

def Open_array_info(filename=''):
    """
    Opening a tiff info, for example size of array, projection and transform matrix.

    Keyword Arguments:
    filename -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file

    """
    try:
        if filename.split('.')[-1] == 'tif':
            f = gdal.Open(r"%s" %filename)
        else:
            f = filename
    except:
            f = filename       
    try:
        geo_out = f.GetGeoTransform()
        proj = f.GetProjection()
        size_X = f.RasterXSize
        size_Y = f.RasterYSize
        f = None
    except:
        print('%s does not exists' %filename)
        
    return(geo_out, proj, size_X, size_Y)

def GDAL_rasterize(input_shp, res, obj):
    
    basefile_shp, ext_shp = os.path.splitext(input_shp)
    output_map = ''.join([basefile_shp, '.tif'])
    name_file = input_shp.split(os.sep)[-1]
    name_File = name_file.split('.')[0]

    # Get environmental variable
    WA_env_paths = os.environ["WA_PATHS"].split(';')
    GDAL_env_path = WA_env_paths[0]
    GDALRASTERIZE_PATH = os.path.join(GDAL_env_path, 'gdal_rasterize.exe') #!!! Vervang deze door gdal.RasterizeLayer
    
    fullCmd = ' '.join(['"%s"' %(GDALRASTERIZE_PATH), '-a %s -tr %s %s -l "%s" "%s" "%s"' %(obj, res, res, name_File, input_shp, output_map)])
    
    Run_command_window(fullCmd)
       
    return(output_map)

def Open_tiff_array(filename='', band=''):
    """
    Opening a tiff array.

    Keyword Arguments:
    filename -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file
    band -- integer
        Defines the band of the tiff that must be opened.
    """
    f = gdal.Open(filename)
    if f is None:
        print('%s does not exists' %filename)
    else:
        if band == '':
            band = 1
        Data = f.GetRasterBand(band).ReadAsArray()
    return(Data)

def Open_nc_info(NC_filename, Var = None):
    """
    Opening a nc info, for example size of array, time (ordinal), projection and transform matrix.

    Keyword Arguments:
    filename -- 'C:/file/to/path/file.nc'
        string that defines the input nc file

    """
    from netCDF4 import Dataset

    fh = Dataset(NC_filename, mode='r')

    if Var is None:
        Var = list(fh.variables.keys())[-1]

    data = fh.variables[Var][:]

    size_Y, size_X = np.int_(data.shape[-2:])
    if len(data.shape) == 3:
        size_Z = np.int_(data.shape[0])
        Time = fh.variables['time'][:]
    else:
        size_Z = 1
        Time = -9999
    lats = fh.variables['latitude'][:]
    lons = fh.variables['longitude'][:]

    Geo6 = fh.variables['latitude'].pixel_size
    Geo2 = fh.variables['longitude'].pixel_size
    Geo4 = np.max(lats) + Geo6/2
    Geo1 = np.min(lons) - Geo2/2

    crso = fh.variables['crs']
    proj = crso.projection
    epsg = Get_epsg(proj, extension = 'GEOGCS')
    geo_out = tuple([Geo1, Geo2, 0, Geo4, 0, -Geo6])
    fh.close()

    return(geo_out, epsg, size_X, size_Y, size_Z, Time)

def Open_nc_array(NC_filename, Var = None, Startdate = '', Enddate = ''):
    """
    Opening a nc array.

    Keyword Arguments:
    filename -- 'C:/file/to/path/file.nc'
        string that defines the input nc file
    Var -- string
        Defines the band name that must be opened.
    Startdate -- "yyyy-mm-dd"
        Defines the startdate (default is from beginning of array)
    Enddate -- "yyyy-mm-dd"
        Defines the enddate (default is from end of array)
    """
    from netCDF4 import Dataset

    fh = Dataset(NC_filename, mode='r')
    if Var == None:
        Var = fh.variables.keys()[-1]

    if Startdate != '':
        Time = fh.variables['time'][:]
        Array_check_start = np.ones(np.shape(Time))
        Date = pd.Timestamp(Startdate)
        Startdate_ord = Date.toordinal()
        Array_check_start[Time >= Startdate_ord] = 0
        Start = np.sum(Array_check_start)
    else:
        Start = 0

    if Enddate != '':
        Time = fh.variables['time'][:]
        Array_check_end = np.zeros(np.shape(Time))
        Date = pd.Timestamp(Enddate)
        Enddate_ord = Date.toordinal()
        Array_check_end[Enddate_ord >= Time] = 1
        End = np.sum(Array_check_end)
    else:
        try:
            Time = fh.variables['time'][:]
            End = len(Time)
        except:
            End = ''

    if (Enddate != '' or Startdate != ''):
        Data = fh.variables[Var][int(Start):int(End), :, :]

    else:
        Data = fh.variables[Var][:]
    fh.close()

    Data = np.array(Data)
    try:
        Data[Data==-9999] = np.nan
    except:
        pass

    return(Data)

def Open_bil_array(bil_filename, band = 1):
    """
    Opening a bil array.

    Keyword Arguments:
    bil_filename -- 'C:/file/to/path/file.bil'
        string that defines the input tiff file or gdal file
    band -- integer
        Defines the band of the tiff that must be opened.
    """    
    gdal.GetDriverByName('EHdr').Register()
    img = gdal.Open(bil_filename)
    Data = img.GetRasterBand(band).ReadAsArray()
    
    return(Data)

def Open_ncs_array(NC_Directory, Var, Startdate, Enddate):
    """
    Opening a nc array.

    Keyword Arguments:
    NC_Directory -- 'C:/file/to/path'
        string that defines the path to all the simulation nc files
    Var -- string
        Defines the band name that must be opened.
    Startdate -- "yyyy-mm-dd"
        Defines the startdate
    Enddate -- "yyyy-mm-dd"
        Defines the enddate
    """

    panda_start = pd.Timestamp(Startdate)
    panda_end = pd.Timestamp(Enddate)

    years = range(int(panda_start.year), int(panda_end.year)+1)
    Data_end = []
    for year in years:

        NC_filename = os.path.join(NC_Directory, "%d.nc" %year)

        if year == years[0]:
            Startdate_now = Startdate
        else:
            Startdate_now = "%d-01-01" %int(year)

        if year == years[-1]:
            Enddate_now = Enddate
        else:
            Enddate_now = "%d-12-31" %int(year)

        Data_now = Open_nc_array(NC_filename, Var, Startdate_now, Enddate_now)

        if year == years[0]:
            Data_end = Data_now
        else:
            Data_end = np.vstack([Data_end,Data_now])

    Data_end = np.array(Data_end)

    return(Data_end)

def Open_nc_dict(input_netcdf, group_name, startdate = '', enddate = ''):
    """
    Opening a nc dictionary.

    Keyword Arguments:
    filename -- 'C:/file/to/path/file.nc'
        string that defines the input nc file
    group_name -- string
        Defines the group name that must be opened.
    Startdate -- "yyyy-mm-dd"
        Defines the startdate (default is from beginning of array)
    Enddate -- "yyyy-mm-dd"
        Defines the enddate (default is from end of array)
    """
    from netCDF4 import Dataset
    import re

    # sort out if the dataset is static or dynamic (written in group_name)
    kind_of_data = group_name.split('_')[-1]

    # if it is dynamic also collect the time parameter
    if kind_of_data == 'dynamic':
        time_dates = Open_nc_array(input_netcdf, Var = 'time')
        Amount_months = len(time_dates)

    # Open the input netcdf and the wanted group name
    in_nc = Dataset(input_netcdf)
    data = in_nc.groups[group_name]

    # Convert the string into a string that can be retransformed into a dictionary
    string_dict = str(data)
    split_dict = str(string_dict.split('\n')[2:-4])
    split_dict = split_dict.replace("'","")
    split_dict = split_dict[1:-1]
    dictionary = dict()
    split_dict_split = re.split(':|,  ',split_dict)

    # Loop over every attribute and add the array
    for i in range(0,len(split_dict_split)):
        number_val = split_dict_split[i]
        if i % 2 == 0:
            Array_text = split_dict_split[i + 1].replace(",","")
            Array_text = Array_text.replace("[","")
            Array_text = Array_text.replace("]","")
            # If the array is dynamic add a 2D array
            if kind_of_data == 'dynamic':
                tot_length = len(np.fromstring(Array_text,sep = ' '))
                dictionary[int(number_val)] = np.fromstring(Array_text,sep = ' ').reshape((int(Amount_months), int(tot_length/Amount_months)))
            # If the array is static add a 1D array
            else:
                dictionary[int(number_val)] = np.fromstring(Array_text,sep = ' ')

    # Clip the dynamic dataset if a start and enddate is defined
    if kind_of_data == 'dynamic':

        if startdate != '':
            Array_check_start = np.ones(np.shape(time_dates))
            Date = pd.Timestamp(startdate)
            Startdate_ord = Date.toordinal()
            Array_check_start[time_dates >= Startdate_ord] = 0
            Start = np.sum(Array_check_start)
        else:
            Start = 0

        if enddate != '':
            Array_check_end = np.zeros(np.shape(time_dates))
            Date = pd.Timestamp(enddate)
            Enddate_ord = Date.toordinal()
            Array_check_end[Enddate_ord >= time_dates] = 1
            End = np.sum(Array_check_end)
        else:
            try:
                time_dates = in_nc.variables['time'][:]
                End = len(time_dates)
            except:
                End = ''

        if Start != 0 or (End != len(time_dates) or ''):

            if End == '':
                End = len(time_dates)

            for key in dictionary.iterkeys():

                Array = dictionary[key][:,:]
                Array_new = Array[int(Start):int(End),:]
                dictionary[key] = Array_new
    in_nc.close()

    return(dictionary)

def Clip_Dataset_GDAL(input_name, output_name, latlim, lonlim):
    """
    Clip the data to the defined extend of the user (latlim, lonlim) by using the gdal_translate executable of gdal.

    Keyword Arguments:
    input_name -- input data, input directory and filename of the tiff file
    output_name -- output data, output filename of the clipped file
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]

    # Get environmental variable
    WA_env_paths = os.environ["WA_PATHS"].split(';')
    GDAL_env_path = WA_env_paths[0]
    GDALTRANSLATE_PATH = os.path.join(GDAL_env_path, 'gdal_translate.exe') #!!! Vervang deze door gdal.Translate

    # find path to the executable    
    fullCmd = ' '.join(["%s" %(GDALTRANSLATE_PATH), '-projwin %s %s %s %s -of GTiff %s %s'  %(lonlim[0], latlim[1], lonlim[1], latlim[0], input_name, output_name)])
    Run_command_window(fullCmd)
    """   
    
    options_list = ['-projwin %s %s %s %s'%(lonlim[0], latlim[1], lonlim[1], latlim[0])]
    options_string = " ".join(options_list)
    gdal.UseExceptions()
    gdal.Translate(output_name, input_name, options = options_string)

    return()

def clip_data(input_file, latlim, lonlim, epsg_clip = 4326, band = 1):
    """
    Clip the data to the defined extend of the user (latlim, lonlim) or to the
    extend of the DEM tile

    Keyword Arguments:
    input_file -- output data, output of the clipped dataset
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    """
    try:
        if input_file.split('.')[-1] == 'tif':
            dest_in = gdal.Open(input_file)
        else:
            dest_in = input_file
    except:
        dest_in = input_file

    # Open Array
    data_in = dest_in.GetRasterBand(band).ReadAsArray()
    Proj_out = Get_epsg(dest_in)
    
    if Proj_out != epsg_clip:
     
        inProj = Proj(init='epsg:%d' %epsg_clip)
        outProj = Proj(init='epsg:%d' %Proj_out)
    
        # Up to here, all  the projection have been defined, as well as a
        # transformation from the from to the to
        lonlim[0], latlim[0] = transform(inProj,outProj,lonlim[0], latlim[0])
        lonlim[1], latlim[1] = transform(inProj,outProj,lonlim[1], latlim[1])       

    # Define the array that must remain
    Geo_in = dest_in.GetGeoTransform()
    Geo_in = list(Geo_in)
    Start_x = np.max([int(np.floor(((lonlim[0]) - Geo_in[0])/ Geo_in[1])),0])
    End_x = np.min([int(np.ceil(((lonlim[1]) - Geo_in[0])/ Geo_in[1])),int(dest_in.RasterXSize)])

    Start_y = np.max([int(np.floor((Geo_in[3] - latlim[1])/ -Geo_in[5])),0])
    End_y = np.min([int(np.ceil(((latlim[0]) - Geo_in[3])/Geo_in[5])), int(dest_in.RasterYSize)])

    #Create new GeoTransform
    Geo_in[0] = Geo_in[0] + Start_x * Geo_in[1]
    Geo_in[3] = Geo_in[3] + Start_y * Geo_in[5]
    Geo_out = tuple(Geo_in)

    data = np.ones([End_y - Start_y, End_x - Start_x])*np.nan

    data = data_in[Start_y:End_y,Start_x:End_x]
    dest_in = None

    return(data, Geo_out, Proj_out)

def clip_data_shp(input_file, shp_filename, IDname = "OBJECTID", number_shp = 1, band = 1):
    """
    Clip the data to the defined extend of the user (latlim, lonlim) or to the
    extend of the DEM tile

    Keyword Arguments:
    input_file -- output data, output of the clipped dataset
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    """
    try:
        if input_file.split('.')[-1] == 'tif':
            dest_in = gdal.Open(input_file)
        else:
            dest_in = input_file
    except:
        dest_in = input_file

    Raster_Basin = Vector_to_Raster(shp_filename, input_file, IDname)
    
    # Create Mask
    MASK = np.where(Raster_Basin == number_shp, 1, np.nan)
    
    # Find IDy
    SumY = np.nansum(MASK, axis = 1)
    IDygood = np.argwhere(SumY>0)
    IDy = [IDygood[0][0], int(IDygood[-1][0]+1)]

    # Find IDX
    SumX = np.nansum(MASK, axis = 0)
    IDxgood = np.argwhere(SumX>0)   
    IDx = [IDxgood[0][0], int(IDxgood[-1][0]+1)]    

    # Get new geotransform
    geo_start = dest_in.GetGeoTransform()
    proj = dest_in.GetProjection()
    
    # Create new geo
    geo_new = tuple([geo_start[0] + geo_start[1] * IDx[0], geo_start[1], 0.0, geo_start[3] + geo_start[5] * IDy[0], 0.0, geo_start[5]])

    # Open Array
    data_in = dest_in.GetRasterBand(band).ReadAsArray()

    # Mask data
    data_temp = data_in * MASK

    # clip data
    data_end = data_temp[IDy[0]:IDy[1],IDx[0]:IDx[1]]

    return(data_end, geo_new, proj)


def reproject_dataset_epsg(dataset, pixel_spacing, epsg_to, method = 2):
    """
    A sample function to reproject and resample a GDAL dataset from within
    Python. The idea here is to reproject from one system to another, as well
    as to change the pixel size. The procedure is slightly long-winded, but
    goes like this:

    1. Set up the two Spatial Reference systems.
    2. Open the original dataset, and get the geotransform
    3. Calculate bounds of new geotransform by projecting the UL corners
    4. Calculate the number of pixels with the new projection & spacing
    5. Create an in-memory raster dataset
    6. Perform the projection

    Keywords arguments:
    dataset -- 'C:/file/to/path/file.tif'
        string that defines the input tiff file
    pixel_spacing -- float
        Defines the pixel size of the output file
    epsg_to -- integer
         The EPSG code of the output dataset
    method -- 1,2,3,4 default = 2
        1 = Nearest Neighbour, 2 = Bilinear, 3 = lanzcos, 4 = average
    """
    import watertools.General.data_conversions as DC

    # 1) Open the dataset
    try:
        g = gdal.Open(dataset)
        if g is None:
            print('input folder does not exist')
    except:
        g = dataset
        
    # Get EPSG code
    epsg_from = Get_epsg(g)

    # Get the Geotransform vector:
    geo_t = g.GetGeoTransform()
    # Vector components:
    # 0- The Upper Left easting coordinate (i.e., horizontal)
    # 1- The E-W pixel spacing
    # 2- The rotation (0 degrees if image is "North Up")
    # 3- The Upper left northing coordinate (i.e., vertical)
    # 4- The rotation (0 degrees)
    # 5- The N-S pixel spacing, negative as it is counted from the UL corner
    x_size = g.RasterXSize  # Raster xsize
    y_size = g.RasterYSize  # Raster ysize

    epsg_to = int(epsg_to)

    # 2) Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
    osng = osr.SpatialReference()
    osng.ImportFromEPSG(epsg_to)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(epsg_from)

    inProj = Proj(init='epsg:%d' %epsg_from)
    outProj = Proj(init='epsg:%d' %epsg_to)

    # Up to here, all  the projection have been defined, as well as a
    # transformation from the from to the to
    ulx, uly = transform(inProj,outProj,geo_t[0], geo_t[3])
    lrx, lry = transform(inProj,outProj,geo_t[0] + geo_t[1] * x_size,
                                        geo_t[3] + geo_t[5] * y_size)

    # See how using 27700 and WGS84 introduces a z-value!
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName('MEM')

    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.
    col = int(np.maximum(int(1), int(np.ceil((lrx - ulx)/pixel_spacing))))
    rows = int(np.maximum(int(1), int(np.ceil((uly - lry)/pixel_spacing))))

    # Re-define lr coordinates based on whole number or rows and columns
    (lrx, lry) = (ulx + col * pixel_spacing, uly -
                  rows * pixel_spacing)

    dest = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)
    dest_nan = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)
    if dest is None:
        print('input folder to large for memory, clip input map')

   # Calculate the new geotransform
    new_geo = (ulx, pixel_spacing, geo_t[2], uly,
               geo_t[4], - pixel_spacing)

    # Set the geotransform
    dest.SetGeoTransform(new_geo)
    dest.SetProjection(osng.ExportToWkt())
    dest_nan.SetGeoTransform(new_geo)
    dest_nan.SetProjection(osng.ExportToWkt())
    
    #
    data_in = g.GetRasterBand(1).ReadAsArray()
    data_in[np.isnan(data_in)] = -9999
    data_nan = np.where(data_in == -9999, 0, 9999)
    g_nan = DC.Save_as_MEM(data_nan, geo_t, epsg_from)

    # Perform the projection/resampling
    if method == 1:
        gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(),gdal.GRA_NearestNeighbour)
    if method == 2:
        gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(),gdal.GRA_Bilinear)
    if method == 3:
        gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Lanczos)
    if method == 4:
        gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Average)
        
    gdal.ReprojectImage(g_nan, dest_nan, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_NearestNeighbour)
    
    data = dest.GetRasterBand(1).ReadAsArray()
    data_nan = dest_nan.GetRasterBand(1).ReadAsArray() 
    data_end = np.where(data_nan != 9999, np.nan, data)
    dest_end = DC.Save_as_MEM(data_end, new_geo, epsg_to)

    return dest_end, ulx, lry, lrx, uly, epsg_to

def reproject_MODIS(input_name, epsg_to):
    '''
    Reproject the merged data file by using gdalwarp. The input projection must be the MODIS projection.
    The output projection can be defined by the user.

    Keywords arguments:
    input_name -- 'C:/file/to/path/file.tif'
        string that defines the input tiff file
    epsg_to -- integer
        The EPSG code of the output dataset
    '''
    
    # Define the output name
    name_out = ''.join(input_name.split(".")[:-1]) + '_reprojected.tif'
   
    src_ds = gdal.Open(input_name)
    
    # Define target SRS
    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromEPSG(int(epsg_to))
    dst_wkt = dst_srs.ExportToWkt()
    
    error_threshold = 0.125  # error threshold --> use same value as in gdalwarp
    resampling = gdal.GRA_NearestNeighbour
    
    # Call AutoCreateWarpedVRT() to fetch default values for target raster dimensions and geotransform
    tmp_ds = gdal.AutoCreateWarpedVRT( src_ds,
                                   None, # src_wkt : left to default value --> will use the one from source
                                   dst_wkt,
                                   resampling,
                                   error_threshold )
    dst_ds = gdal.GetDriverByName('GTiff').CreateCopy(name_out, tmp_ds)
    dst_ds = None 

    return(name_out)

def reproject_shapefile(input_shp, epsg_to):

    input_shp_name_base, input_shp_name_ext= os.path.splitext(input_shp)
    output_shp_name = ''.join([input_shp_name_base, "_EPSG%s" %epsg_to, input_shp_name_ext])
    
    shpf = shapefile.Reader(input_shp)
    wgs_shp = shapefile.Writer(output_shp_name, shapeType=5)
    
    fields = shpf.fields
    
    for name in fields:
        if type(name) == "tuple":
            continue
        else:
            args = name
            wgs_shp.field(*args)
    
    records = shpf.records()
    
    for row in records:
        args = row
        wgs_shp.record(*args)
    
    epsg_from = Get_epsg(input_shp, 'shp')
    
    input_projection = Proj(init="epsg:%s" %epsg_from)
    output_projection = Proj(init="epsg:%s" %epsg_to)
    
    geom = shpf.shapes()
    
    for feature in geom:
        # if there is only one part
        if len(feature.parts) == 1:
            # create empty list to store all the coordinates
            poly_list = []
            # get each coord that makes up the polygon
            for coords in feature.points:
                x, y = coords[0], coords[1]
                # tranform the coord
                new_x, new_y = transform(input_projection, output_projection, x, y)
                # put the coord into a list structure
                poly_coord = [new_x, new_y]
                # append the coords to the polygon list
                poly_list.append(poly_coord)
            # add the geometry to the shapefile.
            wgs_shp.poly([poly_list])
            
        else:
            # append the total amount of points to the end of the parts list
            feature.parts.append(len(feature.points))
            # enpty list to store all the parts that make up the complete feature
            poly_list = []
            # keep track of the part being added
            parts_counter = 0
    
            # while the parts_counter is less than the amount of parts
            while parts_counter < len(feature.parts) - 1:
                # keep track of the amount of points added to the feature
                coord_count = feature.parts[parts_counter]
                # number of points in each part
                no_of_points = abs(feature.parts[parts_counter] - feature.parts[parts_counter + 1])
                # create list to hold individual parts - these get added to poly_list[]
                part_list = []
                # cut off point for each part
                end_point = coord_count + no_of_points
    
                # loop through each part
                while coord_count < end_point:
                    for coords in feature.points[coord_count:end_point]:
                        x, y = coords[0], coords[1]
                        # tranform the coord
                        new_x, new_y = transform(input_projection, output_projection, x, y)
                        # put the coord into a list structure
                        poly_coord = [float(new_x), float(new_y)]
                        # append the coords to the part list
                        part_list.append(poly_coord)
                        coord_count = coord_count + 1
            # append the part to the poly_list
            poly_list.append(part_list)
            parts_counter = parts_counter + 1
            # add the geometry to to new file
            wgs_shp.poly([poly_list])   
    
    
    wgs_shp.close()

    prj = open(''.join([os.path.splitext(output_shp_name)[0],".prj"]), "w")
    epsg = getWKT_PRJ("%s" %epsg_to)
    prj.write(epsg)
    prj.close()
    
    return(output_shp_name)

def reproject_dataset_example(dataset, dataset_example, method=1):
    """
    A sample function to reproject and resample a GDAL dataset from within
    Python. The user can define the wanted projection and shape by defining an example dataset.

    Keywords arguments:
    dataset -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file
    dataset_example -- 'C:/file/to/path/file.tif' or a gdal file (gdal.Open(filename))
        string that defines the input tiff file or gdal file
    method -- 1,2,3,4 default = 1
        1 = Nearest Neighbour, 2 = Bilinear, 3 = lanzcos, 4 = average
    """
    # open dataset that must be transformed
    try:
        if (os.path.splitext(dataset)[-1] == '.tif' or os.path.splitext(dataset)[-1] == '.TIF'):
            g = gdal.Open(dataset)
        else:
            g = dataset
    except:
            g = dataset
    epsg_from = Get_epsg(g)

    #exceptions
    if epsg_from == 9001:
        epsg_from = 5070

    # open dataset that is used for transforming the dataset
    try:
        if (os.path.splitext(dataset_example)[-1] == '.tif' or os.path.splitext(dataset_example)[-1] == '.TIF'):
            gland = gdal.Open(dataset_example)
            epsg_to = Get_epsg(gland)
        elif os.path.splitext(dataset_example)[-1] == '.nc':
            import watertools.General.data_conversions as DC
            geo_out, epsg_to, size_X, size_Y, size_Z, Time = Open_nc_info(dataset_example)
            data = np.zeros([size_Y, size_X])
            gland = DC.Save_as_MEM(data, geo_out, str(epsg_to))
        else:
            gland = dataset_example
            epsg_to = Get_epsg(gland)
    except:
            gland = dataset_example
            epsg_to = Get_epsg(gland)

    # Set the EPSG codes
    osng = osr.SpatialReference()
    osng.ImportFromEPSG(epsg_to)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(epsg_from)

    # Get shape and geo transform from example
    geo_land = gland.GetGeoTransform()
    col=gland.RasterXSize
    rows=gland.RasterYSize

    # Create new raster
    mem_drv = gdal.GetDriverByName('MEM')
    dest1 = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)
    dest1.SetGeoTransform(geo_land)
    dest1.SetProjection(osng.ExportToWkt())

    # Perform the projection/resampling
    if method == 1:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_NearestNeighbour)
    if method == 2:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Bilinear)
    if method == 3:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Lanczos)
    if method == 4:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Average)
    if method == 5:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Cubic)
    if method == 6:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_CubicSpline)
    if method == 7:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Mode)
    if method == 8:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Max)
    if method == 9:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Min)
    if method == 10:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Med)
    if method == 11:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Q1)
    if method == 12:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Q3)
    if method == 13:
        gdal.ReprojectImage(g, dest1, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Sum)
        
    return(dest1)

def resize_array_example(Array_in, Array_example, method=1):
    """
    This function resizes an array so it has the same size as an example array
    The extend of the array must be the same

    Keyword arguments:
    Array_in -- []
        Array: 2D or 3D array
    Array_example -- []
        Array: 2D or 3D array
    method: -- 1 ... 5
        int: Resampling method
    """

    # Create old raster
    Array_out_shape = np.int_(Array_in.shape)
    Array_out_shape[-1] = Array_example.shape[-1]
    Array_out_shape[-2] = Array_example.shape[-2]

    if method == 1:
        interpolation_method='nearest'
        interpolation_number = 0
    if method == 2:
        interpolation_method='bicubic'
        interpolation_number = 3
    if method == 3:
        interpolation_method='bilinear'
        interpolation_number = 1   
    if method == 4:
        interpolation_method='cubic'     
    if method == 5:
        interpolation_method='lanczos'

    if len(Array_out_shape) == 3:
        Array_out = np.zeros(Array_out_shape)

        for i in range(0, Array_out_shape[0]):
            Array_in_slice = Array_in[i,:,:]
            size=tuple(Array_out_shape[1:])

            if sys.version_info[0] == 2:
                import scipy.misc as misc
                Array_out_slice= misc.imresize(np.float_(Array_in_slice), size, interp=interpolation_method, mode='F')
            if sys.version_info[0] == 3:
                import skimage.transform as transform
                Array_out_slice= transform.resize(np.float_(Array_in_slice), size, order=interpolation_number)
                
            Array_out[i,:,:] = Array_out_slice

    elif len(Array_out_shape) == 2:

        size=tuple(Array_out_shape)
        if sys.version_info[0] == 2:
            import scipy.misc as misc
            Array_out= misc.imresize(np.float_(Array_in), size, interp=interpolation_method, mode='F')
        if sys.version_info[0] == 3:
            import skimage.transform as transform
            Array_out= transform.resize(np.float_(Array_in), size, order=interpolation_number)
      
    else:
        print('only 2D or 3D dimensions are supported')

    return(Array_out)

def Get_epsg(g, extension = 'tiff'):
    """
    This function reads the projection of a GEOGCS file or tiff file

    Keyword arguments:
    g -- string
        Filename to the file that must be read
    extension -- tiff or GEOGCS
        Define the extension of the dataset (default is tiff)
    """
    try:
        if extension == 'tiff':
            # Get info of the dataset that is used for transforming
            try:
                dest = gdal.Open(g)
            except:
                dest = g
            g_proj = dest.GetProjection()
            Projection=g_proj.split('EPSG","')
            epsg_to=int((str(Projection[-1]).split(']')[0])[0:-1])
                
        if extension == 'GEOGCS':
            Projection = g
            epsg_to=int((str(Projection).split('"EPSG","')[-1].split('"')[0:-1])[0])
        if extension == 'shp':
            data = fiona.open(g)
            epsg_to = int(data.crs['init'].split(":")[-1])
            
            '''
            projection_file = ''.join([os.path.splitext(g)[0],'.prj'])
            print(projection_file)
            Projection = open(projection_file, 'r').read()
            print(Projection)            
            srs = osr.SpatialReference()
            srs.SetFromUserInput(Projection)
            epsg_to = srs.GetAttrValue("AUTHORITY",1)
            print(srs)
            if epsg_to == None:
                try:
                    epsg_str = srs.GetAttrValue("PROJCS", 0)        
                    zone = epsg_str.split("_")[-1][0:-1]
                    NorS = str(epsg_str.split("_")[-1][-1])
                    if NorS == "N":
                        SN = 6
                    if NorS == "S":
                        SN = 7   
                    epsg_to = int("32%s%02s" %(SN, zone))
                    
                except:       
                    epsg_to=4326   
            '''                    
    except:
        epsg_to=4326
        #print 'Was not able to get the projection, so WGS84 is assumed'
        
    return(epsg_to)

def gap_filling(dataset, NoDataValue, method = 1):
    """
    This function fills the no data gaps in a numpy array

    Keyword arguments:
    dataset -- 'C:/'  path to the source data (dataset that must be filled)
    NoDataValue -- Value that must be filled
    """
    import watertools.General.data_conversions as DC

    try:
        if dataset.split('.')[-1] == 'tif':
            # Open the numpy array
            data = Open_tiff_array(dataset)
            Save_as_tiff = 1
        else:
            data = dataset
            Save_as_tiff = 0
    except:
        data = dataset
        Save_as_tiff = 0

    # fill the no data values
    if NoDataValue == np.nan:
        mask = ~(np.isnan(data))
    else:
        mask = ~(data==NoDataValue)
    xx, yy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    xym = np.vstack( (np.ravel(xx[mask]), np.ravel(yy[mask])) ).T
    data0 = np.ravel( data[:,:][mask] )

    if method == 1:
        interp0 = scipy.interpolate.NearestNDInterpolator( xym, data0 )
        data_end = interp0(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )

    if method == 2:
        interp0 = scipy.interpolate.LinearNDInterpolator( xym, data0 )
        data_end = interp0(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )

    if Save_as_tiff == 1:
        EndProduct=dataset[:-4] + '_GF.tif'

        # collect the geoinformation
        geo_out, proj, size_X, size_Y = Open_array_info(dataset)

        # Save the filled array as geotiff
        DC.Save_as_tiff(name=EndProduct, data=data_end, geo=geo_out, projection=proj)

    else:
        EndProduct = data_end

    return (EndProduct)

def Get3Darray_time_series_monthly(Data_Path, Startdate, Enddate, Example_data = None):
    """
    This function creates a datacube

    Keyword arguments:
    Data_Path -- 'product/monthly'
        str: Path to the dataset
    Startdate -- 'YYYY-mm-dd'
        str: startdate of the 3D array
    Enddate -- 'YYYY-mm-dd'
        str: enddate of the 3D array
    Example_data: -- 'C:/....../.tif'
        str: Path to an example tiff file (all arrays will be reprojected to this example)
    """

    # Get a list of dates that needs to be reprojected
    Dates = pd.date_range(Startdate, Enddate, freq = 'MS')

    # Change Working directory
    os.chdir(Data_Path)
    i = 0

    # Loop over the months
    for Date in Dates:

        # Create the end monthly file name
        End_tiff_file_name = 'monthly_%d.%02d.01.tif' %(Date.year, Date.month)

        # Search for this file in directory
        file_name = glob.glob('*%s' %End_tiff_file_name)

        # Select the first file that is found
        file_name_path = os.path.join(Data_Path, file_name[0])

        # Check if an example file is selected
        if Example_data is not None:

            # If it is the first day set the example gland file
            if Date == Dates[0]:

                # Check the format to read general info

                # if Tiff
                if os.path.splitext(Example_data)[-1] == '.tif':
                    geo_out, proj, size_X, size_Y = Open_array_info(Example_data)
                    dataTot=np.zeros([len(Dates),size_Y,size_X])

                # if netCDF
                if os.path.splitext(Example_data)[-1] == '.nc':
                    geo_out, projection, size_X, size_Y, size_Z, Time = Open_nc_info(Example_data)
                    dataTot=np.zeros([len(Dates),size_Y,size_X])

                    # Create memory file for reprojection
                    data = Open_nc_array(Example_data, "Landuse")
                    driver = gdal.GetDriverByName("MEM")
                    gland = driver.Create('', int(size_X), int(size_Y), 1,
                                           gdal.GDT_Float32)
                    srse = osr.SpatialReference()
                    if projection == '' or projection == 4326:
                        srse.SetWellKnownGeogCS("WGS84")
                    else:
                        srse.SetWellKnownGeogCS(projection)
                    gland.SetProjection(srse.ExportToWkt())
                    gland.GetRasterBand(1).SetNoDataValue(-9999)
                    gland.SetGeoTransform(geo_out)
                    gland.GetRasterBand(1).WriteArray(data)

                # use the input parameter as it is already an example file
                else:
                    gland = Example_data

            # reproject dataset
            dest = reproject_dataset_example(file_name_path, gland, method = 4)
            Array_one_date = dest.GetRasterBand(1).ReadAsArray()

        # if there is no example dataset defined
        else:

            # Get the properties from the first file
            if Date is Dates[0]:
                    geo_out, proj, size_X, size_Y = Open_array_info(file_name_path)
                    dataTot=np.zeros([len(Dates),size_Y,size_X])
            Array_one_date = Open_tiff_array(file_name_path)

        # Create the 3D array
        dataTot[i,:,:] = Array_one_date
        i += 1

    return(dataTot)

def Vector_to_Raster(shapefile_name, filename, Attribute_name):
    """
    This function creates a raster of a shp file

    Keyword arguments:
    shapefile_name -- 'C:/....../.shp'
        str: Path from the shape file
    reference_raster_data_name -- destination file as example file

    """
    from osgeo import gdal, ogr

    try:
        if filename.split('.')[-1] == 'tif':
            dest_ex = gdal.Open(r"%s" %filename)
        else:
            dest_ex = filename
    except:
            dest_ex = filename      

    proj = Get_epsg(dest_ex)
    geo = dest_ex.GetGeoTransform()
    size_X = dest_ex.RasterXSize
    size_Y = dest_ex.RasterYSize
    
    x_min = geo[0]
    x_max = geo[0] + size_X * geo[1]
    y_min = geo[3] + size_Y * geo[5]
    y_max = geo[3]
    pixel_size_x = geo[1]
    pixel_size_y = -geo[5]
    
    # Open the data source and read in the extent
    source_ds = ogr.Open(shapefile_name)
    source_layer = source_ds.GetLayer()

    # Create the destination data source
    x_res = int(round((x_max - x_min) / pixel_size_x))
    y_res = int(round((y_max - y_min) / pixel_size_y))

    # Create tiff file
    target_ds = gdal.GetDriverByName('MEM').Create('', x_res, y_res, 1, gdal.GDT_Float32)
    target_ds.SetGeoTransform(geo)
    srse = osr.SpatialReference()
    srse.ImportFromEPSG(int(proj))
    target_ds.SetProjection(srse.ExportToWkt())
    band = target_ds.GetRasterBand(1)
    target_ds.GetRasterBand(1).SetNoDataValue(-9999)
    band.Fill(-9999)

    # Rasterize the shape and save it as band in tiff file
    gdal.RasterizeLayer(target_ds, [1], source_layer, options=["ATTRIBUTE=%s" %Attribute_name])

    # Open array
    Raster_Basin = target_ds.GetRasterBand(1).ReadAsArray()
    target_ds = None
    
    return(Raster_Basin)

def Moving_average(dataset, Moving_front, Moving_back):
    """
    This function applies the moving averages over a 3D matrix called dataset.

    Keyword Arguments:
    dataset -- 3D matrix [time, ysize, xsize]
    Moving_front -- Amount of time steps that must be considered in the front of the current month
    Moving_back -- Amount of time steps that must be considered in the back of the current month
    """

    dataset_out = np.zeros((int(np.shape(dataset)[0]) - Moving_back - Moving_front, int(np.shape(dataset)[1]), int(np.shape(dataset)[2])))

    for i in range(Moving_back, (int(np.shape(dataset)[0]) - Moving_front)):
        dataset_out[i - Moving_back,:,:] = np.nanmean(dataset[i - Moving_back : i + 1 + Moving_front, :,:], 0)

    return(dataset_out)

def Get_ordinal(Startdate, Enddate, freq = 'MS'):
    """
    This function creates an array with ordinal time.

    Keyword Arguments:
    Startdate -- Startdate of the ordinal time
    Enddate -- Enddate of the ordinal time
    freq -- Time frequencies between start and enddate
    """

    import datetime
    Dates = pd.date_range(Startdate, Enddate, freq = freq)
    i = 0
    ordinal = np.zeros([len(Dates)])
    for date in Dates:

        p = datetime.date(date.year, date.month, date.day).toordinal()
        ordinal[i]=p
        i += 1

    return(ordinal)

def Create_Buffer(Data_In, Buffer_area = 2):

   '''
   This function creates a 3D array which is used to apply the moving window
   '''

   # Buffer_area = 2 # A block of 2 times Buffer_area + 1 will be 1 if there is the pixel in the middle is 1
   Data_Out=np.empty((len(Data_In),len(Data_In[1])))
   Data_Out[:,:] = Data_In
   for ypixel in range(0,Buffer_area + 1):

        for xpixel in range(1,Buffer_area + 1):

           if ypixel==0:
                for xpixel in range(1,Buffer_area + 1):
                    Data_Out[:,0:-xpixel] += Data_In[:,xpixel:]
                    Data_Out[:,xpixel:] += Data_In[:,:-xpixel]

                for ypixel in range(1,Buffer_area + 1):

                    Data_Out[ypixel:,:] += Data_In[:-ypixel,:]
                    Data_Out[0:-ypixel,:] += Data_In[ypixel:,:]

           else:
               Data_Out[0:-xpixel,ypixel:] += Data_In[xpixel:,:-ypixel]
               Data_Out[xpixel:,ypixel:] += Data_In[:-xpixel,:-ypixel]
               Data_Out[0:-xpixel,0:-ypixel] += Data_In[xpixel:,ypixel:]
               Data_Out[xpixel:,0:-ypixel] += Data_In[:-xpixel,ypixel:]

   Data_Out[Data_Out>0.1] = 1
   Data_Out[Data_Out<=0.1] = 0

   return(Data_Out)
   
def getWKT_PRJ(epsg_code):   

    if epsg_code == str(900913) or epsg_code == str(3857):
        output = 'PROJCS["WGS_1984_Web_Mercator_Auxiliary_Sphere",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Mercator_Auxiliary_Sphere"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["Standard_Parallel_1",0.0],PARAMETER["Auxiliary_Sphere_Type",0.0],UNIT["Meter",1.0]]'
    
    else:
        if sys.version_info[0] == 3:
            wkt = requests.get(
                "http://spatialreference.org/ref/epsg/{0}/prettywkt/"
                .format(epsg_code))
            
            if wkt.status_code == 502:
                print("http://spatialreference.org is offline, osr is used to get"\
                         "projection")

                srs = osr.SpatialReference()
                srs.ImportFromEPSG(int(epsg_code))
                output = srs.ExportToWkt()
                
            else:
                input_str = wkt.content.decode('utf-8')
                remove_space = input_str.replace(" ","")     
                output = remove_space.replace("\n", "")                    
                        
    
        if sys.version_info[0] == 2:
            import urllib2
            wkt = urllib2.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg_code))
            remove_spaces = wkt.read().replace(" ","")
            output = remove_spaces.replace("\n", "")   

    return output   


def Update_log_website(url, postfields):
    
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.POSTFIELDS, postfields)
    c.perform()    
    c.close()
    
    return()   