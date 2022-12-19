# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 13:14:04 2022

@author: timhe
"""
import sys
import re
import datetime
import os
import requests
import pycurl
import numpy as np
import pandas as pd
import urllib
from osgeo import gdal
from netCDF4 import Dataset
import shutil

import watertools
import watertools.General.data_conversions as DC
import watertools.General.raster_conversions as RC

# input_folder = r'H:/Project_WAKAS_test/Data_Out/Input_Data'

# startdate = '2022-11-25'
# enddate = '2022-12-27'
# latlim = [29.8, 31.3]
# lonlim = [71.5, 74]

def DownloadData(input_folder, startdate, enddate, latlim, lonlim, Waitbar):
    
    # Estimate GMT Offset
    lon = np.nanmean(lonlim)
    GMT_Offset = round(lon * 24 / 360)
    
    Dates = pd.date_range(startdate, enddate, freq = "D")
    
    input_folder_VIIRS = os.path.join(input_folder, "VIIRS", "LST", "Daily")
    input_folder_VIIRS_RAW = os.path.join(input_folder, "VIIRS", "LST", "Daily", "RAW")
    
    if not os.path.exists(input_folder_VIIRS_RAW):
        os.makedirs(input_folder_VIIRS_RAW)
        
    # Create Waitbar
    total_amount = len(Dates)
    if Waitbar == 1:
        import watertools.Functions.Random.WaitbarConsole as WaitbarConsole
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    for Date in Dates:
        
        date_str = datetime.datetime.strftime(Date, "%Y-%m-%d")
        
        # Search for content
        firstPage()
        content_viirs = Seach_For_VIIRS(date_str, date_str, latlim, lonlim, GMT_Offset)
        
        # Find filename GITCO and SVI
        print("Find filenames GITCO and SVI")
        filenames_SVI, filenames_GITCO = Get_Filenames_VIIRS(content_viirs)
       
        for filename_svi in filenames_SVI:
                
            print("Download RAW data for %s" %filename_svi)
            
            date_datetime =  datetime.datetime.strptime(filename_svi.split("_")[2][1:], "%Y%m%d")         
            time_str_start = filename_svi.split("_")[3][1:]
            time_str_end = filename_svi.split("_")[4][1:]
            DOY = date_datetime.strftime("%j")
            
            url_viaes_start = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/NPP_VIAES_L1/%s/%s/" %(date_datetime.year, DOY)
            url_vnp_start = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/VNP03IMG/%s/%s/" %(date_datetime.year, DOY)
        
            # Find right files
            r = requests.get("%s.csv" %url_viaes_start)
            files_viaes = str(r.content)
            files_viaes_hdf = [files_viaes[m.start():m.start()+48] for m in re.finditer("NPP_VIAES_L1.A", files_viaes)]
            files_viaes_hdf_ok = [k for k in files_viaes_hdf if int(k.split(".")[2])>=int(time_str_start[0:4]) - 4 and int(k.split(".")[2])<=int(time_str_end[0:4]) + 4]
            
            r = requests.get("%s.csv" %url_vnp_start)
            files_vnp = str(r.content)                
            files_vnp_hdf = [files_vnp[m.start():m.start()+43] for m in re.finditer("VNP03IMG.A", files_vnp)]
            files_vnp_hdf_ok = [k for k in files_vnp_hdf if int(k.split(".")[2])>=int(time_str_start[0:4]) - 4 and int(k.split(".")[2])<=int(time_str_end[0:4]) + 4]
            
            if len(files_viaes_hdf_ok) > 0 and len(files_vnp_hdf_ok) > 0:
                
                for idx, file_viaes_hdf_ok in enumerate(files_viaes_hdf_ok):
                    
                    file_vnp_hdf_ok = files_vnp_hdf_ok[idx]
                    
                    print("OK VIEAS file: %s" %file_viaes_hdf_ok)
                    print("OK VNP file: %s" %file_vnp_hdf_ok)         
                    
                    # Check if time is same:
                    if int(file_viaes_hdf_ok.split(".")[2][0:4]) == int(file_vnp_hdf_ok.split(".")[2][0:4]): 
                        
                        hour = int(file_viaes_hdf_ok.split(".")[2][0:2])
                        minutes = int(file_viaes_hdf_ok.split(".")[2][2:4]) 
                        
                        filename_out = os.path.join(input_folder_VIIRS, "NPP_VIAES_L1_%s_H%02dM%02d.tif" %(date_datetime.strftime("%Y%m%d"), hour, minutes))
                        
                        if not os.path.exists(filename_out):
                            
                            file_name_viaes = os.path.join(input_folder_VIIRS_RAW, file_viaes_hdf_ok)                
                            file_name_vnp = os.path.join(input_folder_VIIRS_RAW, file_vnp_hdf_ok)
                                                      
                            try:
                                
                                 print("Download VNP03IMG from http %s" %file_vnp_hdf_ok)
                                 # Download VNP03IMG
                                 url_vnp = url_vnp_start + file_vnp_hdf_ok
                                 Perform_Download_VIIRS(url_vnp, file_name_vnp)
                                 
                                 # Check if data is required
                                 fh2 = Dataset(file_name_vnp)   
                                 
                                 # Open geo
                                 lat = fh2['geolocation_data']['latitude'][:,:].data
                                 lon = fh2['geolocation_data']['longitude'][:,:].data
                                 
                                 # Check if data points are in AOI
                                 lat_data = lat.flatten()
                                 lon_data = lon.flatten()
                                
                                 # make search area a bit larger
                                 latlim_OK = [latlim[0] - 0.1, latlim[1] + 0.1]
                                 lonlim_OK = [lonlim[0] - 0.1, lonlim[1] + 0.1]    
                
                                 boolean = np.where(np.logical_and.reduce([lat_data >= latlim_OK[0], latlim_OK[1] >= lat_data, lon_data >=lonlim_OK[0], lonlim_OK[1] >= lon_data]), 1, 0)
                                
                                 lat_data_ok = lat_data[boolean==1]
                                 lon_data_ok = lon_data[boolean==1]    
                                 fh2.close()
                                 
                                 if len(lat_data_ok) > 0:
                                     
                                     print("Download NPP_VIEAS_L1 from http %s" %file_viaes_hdf_ok)
                                     
                                     attemps = 0
                                     success = 0
                                     while attemps<4 and success==0:
                                         
                                         # Download NPP_VIEAS_L1
                                         url_viaes = url_viaes_start + file_viaes_hdf_ok
                                         Perform_Download_VIIRS(url_viaes, file_name_viaes)
                                         attemps +=1
                                         stats = os.stat(file_name_viaes)
                                         
                                         if stats.st_size > 10000:
                                             success = 1

                                     print("Combine VNP03IMG and NPP_VIEAS_L1")
                                     
                                    
                                     # Open data  
                                     fh = gdal.Open(file_name_viaes)
                                     sdsdict = fh.GetSubDatasets()
                                     fh = gdal.Open(sdsdict[13][0])
                                    
                                     # Open LST
                                     LST = fh.ReadAsArray()
                                     LST = np.float_(LST)
                                     LST[LST>60000] = np.nan
                                     LST = 0.003510003444 * LST + 150
                                     del fh
        
                                     print("Select the right points for an AOI that is a bit larger")
                                     # Flatten data   
                                     z_data = LST.flatten()
                                     z_data_ok = z_data[boolean==1]   
                                                           
                                        
                                     dfn = pd.DataFrame({"x":lon_data_ok, "y":lat_data_ok, "value":z_data_ok})
                                     dfn2 = dfn.dropna(subset=['value'])
                                    
                                     os.chdir(input_folder_VIIRS_RAW)
                                     csv_out = os.path.join(input_folder_VIIRS_RAW, "dfn.csv")
                                     dfn2.to_csv(csv_out, index = False, sep = ",")
                                    
                                     vrt_out = os.path.join(input_folder_VIIRS_RAW, "dfn.vrt") 
                                     tif_out = os.path.join(input_folder_VIIRS_RAW, "dfn.tif") 
                                    
                                     print("Create virtual layer")
                                     f = open(vrt_out, "w")
                                     f.write("<OGRVRTDataSource>\n \
                                         <OGRVRTLayer name=\"dfn\">\n \
                                             <SrcDataSource>dfn.csv</SrcDataSource>\n \
                                             <GeometryType>wkbPoint</GeometryType>\n \
                                             <GeometryField encoding=\"PointFromColumns\" x=\"x\" y=\"y\" z=\"value\"/>\n \
                                         </OGRVRTLayer>\n \
                                     </OGRVRTDataSource>")
                                    
                                     f.close()
                                    
                                     print("Rasterize data")
                                     dest = gdal.Rasterize(tif_out, vrt_out, outputSRS = "EPSG:4326", xRes = 0.00375, yRes = -0.00375, attribute = "value", noData = np.nan)
                                    
                                     LST = dest.GetRasterBand(1).ReadAsArray()
                                     LST[LST==0] = np.nan
                                     
                                     BUFFER = np.where(~np.isnan(LST), 1, 0)
                                     BUFFER_ONE = RC.Create_Buffer(BUFFER, 1)
                                     
                                     geo = dest.GetGeoTransform()
                                     proj = dest.GetProjection()
                                     dest = None
                                                
                                     print("Apply gapfilling")
                                     LST[np.logical_and(np.isnan(LST), BUFFER_ONE==1)] = -9999
                                     #LST_GF = RC.gap_filling(LST, -9999, method = 2)
                                     #LST_GF[np.isnan(LST_GF)] = -9999
                                     LST_GF = RC.gap_filling(LST, -9999, method = 1)
                                     
                                     print("Create %s" %filename_out)
                                     DC.Save_as_tiff(filename_out, LST_GF, geo, proj)   
            
                            except:
                                
                                print("Error creating VIIRS file for %s" %file_viaes_hdf_ok)
        
        # Adjust waitbar
        if Waitbar == 1:
            import watertools.Functions.Random.WaitbarConsole as WaitbarConsole
            amount += 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
        
    try:
        shutil.rmtree(input_folder_VIIRS_RAW)        
    except:
        print("Was not able to delete %s" %input_folder_VIIRS_RAW)
    
    return()


    
def Perform_Download_VIIRS(url, file_name_out):
    
    # url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/VNP03IMG/2022/334/VNP03IMG.A2022334.0748.002.2022334152203.nc'
    # file_name_out = 'H:/Project_WAKAS_test/Data_Out/Input_Data\\VIIRS\\LST\\Daily\\RAW\\VNP03IMG.A2022334.0748.002.2022334152203.nc'
    
    VIIRS_BEARER, passw = watertools.Functions.Random.Get_Username_PWD.GET('VIIRS_BEARER')
    head = {"Authorization": "Bearer %s" %VIIRS_BEARER}
    
    r = requests.get(url, headers = head, stream=True, timeout = 2000)
    
    with open(file_name_out, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
                
    return()

def firstPage(site = "avl"):
    firstPage = 'https://www.%s.class.noaa.gov/saa/products/classlogin' %site
    if os.path.exists('cookie.txt'):
        os.remove('cookie.txt')

    try:
        c = pycurl.Curl()
        c.setopt(c.URL, firstPage)
        c.setopt(pycurl.COOKIEJAR, 'cookie.txt')
        c.setopt(pycurl.COOKIEFILE, 'cookie.txt')
        c.setopt(pycurl.TIMEOUT, 1200)        
        s = c.perform_rs()
        c.close()
    except Exception as E:
        print("ERROR: %s" %E) 
        print('Unable to fetch firstPage')


def Seach_For_VIIRS(date_start, date_end, latlim, lonlim, GMT_Offset):
    
    print("Search on site for filenames")    
    putInCartPage = 'https://www.avl.class.noaa.gov/saa/products/psearchVIIRS_SDR'
    params = {}
    params['search_opt'] = 'SC'
    params['gid_pattern'] = '^(NPP|J01)\d{12}$'
    params['orb_pattern'] = '^(\d{1,})$'
    params['dsname_pattern'] = '^((\w{5})_(NPP|J0[12])|\w{5}-\w{5}_(NPP|J0[12]))_D20\d\d(0[1-9]|1[012])([012][0-9]|3[01]).*$'
    params['nlat'] = latlim[1]
    params['wlon'] = lonlim[0]
    params['elon'] = lonlim[1]
    params['slat'] = latlim[0]
    params['minDiff'] = '0.0'
    params['data_start'] = '2008-04-30'
    params['data_end'] = date_end
    params['max_days_val'] = '366'
    params['start_date'] = date_start
    params['start_time'] = '%s:00:00' %int(np.floor(9-GMT_Offset))
    params['end_date'] = date_end
    params['end_time'] = '%s:00:00' %int(np.ceil(17-GMT_Offset))
    params['between_through'] = 'T'
    params['max_sum_hits'] = '2000'
    params['lrg_max_sum_hits'] = '4000'
    params['brk_srch_hrs_qs'] = '6'
    params['bulk_order'] = 'N'
    params['limit_search'] = 'Y'
    params['max_lat_range'] = '180'
    params['max_lon_range'] = '360'
    params['datatype_family'] = 'VIIRS_SDR'

    params['Datatype'] = ['VIIRSI5SDR', 'VIIRSITGEO']
    params['Satellite'] = 'NPP'
    params['Node'] = 'A'

    print("VIIRS search params: %s" %params)

    # Single valued items
    safeCharacters = '()*'
    paramList = urllib.parse.urlencode([(h,v) for (h,v) in params.items() if not isinstance(v,list)], safe=safeCharacters)
    # Multivalued items
    paramList += "&"
    paramList += "&".join(["&".join([h+'='+urllib.parse.quote(v, safe=safeCharacters) for v in v2]) for (h,v2) in params.items() if isinstance(v2, list)])

    try:
        c = pycurl.Curl()
        c.setopt(c.POST, 1)
        c.setopt(c.URL, putInCartPage)
        c.setopt(pycurl.FOLLOWLOCATION, 1)
        c.setopt(pycurl.COOKIEJAR, 'cookie.txt')
        c.setopt(pycurl.COOKIEFILE, 'cookie.txt')
        c.setopt(pycurl.POSTFIELDS, paramList)
        c.setopt(pycurl.TIMEOUT, 1200)
        s = c.perform_rs()
    except:
        print("ERROR:", sys.exc_info()[0], sys.exc_info()[1])
        s = ""
    print("VIIRS search finished")
    return(s)
    
def Get_Filenames_VIIRS(content_viirs):
    
    #print(content_viirs)
    filenames_SVI = []
    SVI_ends = content_viirs.split("SVI05_npp_d")
    
    if len(SVI_ends)>1:
        for SVI_end in SVI_ends[1:]:
            SVI_end = SVI_end.split(".h5")
            filenames_SVI.append("SVI05_npp_d" + SVI_end[0].replace("<wbr/>", "") + ".h5")

    filenames_GITCO = []        
    GITCO_ends = content_viirs.split("GITCO_npp_d")
    
    if len(GITCO_ends)>1:
        for GITCO_end in GITCO_ends[1:]:
            GITCO_end = GITCO_end.split(".h5")
            filenames_GITCO.append("GITCO_npp_d" + GITCO_end[0].replace("<wbr/>", "") + ".h5")
        
    return(filenames_SVI, filenames_GITCO )

