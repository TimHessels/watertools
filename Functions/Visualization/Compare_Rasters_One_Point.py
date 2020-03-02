# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 11:34:43 2018

@author: tih
"""
import numpy as np
import pandas as pd
import glob
import os
import gdal
import osr
from pyproj import Proj, transform

#input_folders = [r"G:\SEBAL_Tadla\Daily_SEBAL\Dataset_out\PreSEBAL_SEBAL_out\NDVI",r"G:\SEBAL_Tadla\Daily_SEBAL\Dataset_out\NDVI_daily",""]
#input_formats = ["*{yyyy}{mm:02d}{dd:02d}.tif", "*{yyyy}_{doy}.tif"]
#Coordinate = [-6.80689, 32.46206]
#Startdate = "2016-10-01"
#Enddate = "2017-09-30"

def Visualize_Graph(input_folders, input_formats, Coordinate, Startdate, Enddate, freq = "D", vmin = None, vmax = None):

    Dates = pd.date_range(Startdate, Enddate, freq = freq)
    
    Datasets_end = np.ones([len(Dates), len(input_folders)+1]) * np.nan
    
    for j in range(0,len(input_folders)):
    
        input_folder = input_folders[j]  
        input_format = input_formats[j]
        Dataset = Get_Dataset_Point(input_folder, input_format, Dates, Coordinate)    
        
        if j == 0:
            
            Datasets_end[:,0:2] = Dataset
            
        else:
            Datasets_end[:,j + 1] = Dataset[:,1]
            
    Datasets_end[Datasets_end==-9999]=np.nan
    
    import matplotlib.pyplot as plt
    
    colors = ['red','blue','green','black','yellow']
    
    for k in range(0,len(input_folders)):
        plt.scatter(Dates, Datasets_end[:,k+1], s =10, c = colors[k])
    plt.ylim(vmin, vmax)
    plt.legend(input_folders, loc='upper center',bbox_to_anchor=(0.5, -0.1))
    plt.show()  
    
    return(plt)



def Get_Dataset_Point(input_folder, input_format, Dates, Coordinate):
      
    Dataset = np.ones([len(Dates), 2]) * np.nan
    i = 0
    
    for Date in Dates:

        Day = Date.day
        Month = Date.month
        Year = Date.year
        Hour = Date.hour
        Minute = Date.minute
        DOY = Date.dayofyear
        Dataset[i,0] = '%d%02d%02d' %(Year, Month, Day)
        os.chdir(input_folder)
        filename = glob.glob(input_format.format(yyyy=Year, mm=Month, dd=Day, doy = DOY, HH = Hour, MM=Minute))
    
        if len(filename) > 0: 
            try:
                filename = filename[0]
                
                if not'yID' in locals():
                    yID, xID = Get_Row_Column_CoordinateWGS(filename, Coordinate)
                
                if (np.isnan(yID) or np.isnan(xID)):
                    Value = np.nan
                else:
                    dest = gdal.Open(filename)
                    Array = dest.GetRasterBand(1).ReadAsArray()
                    Value = Array[yID, xID]
                
                Dataset[i,1] = Value
            except:   
                Dataset[i,1] = np.nan
        i += 1

    del yID, xID
    
    return(Dataset)

    
def Get_Row_Column_CoordinateWGS(filename, Coordinate):
             
    dest = gdal.Open(filename)
    geo_out = dest.GetGeoTransform()
    size_x = dest.RasterXSize
    size_y = dest.RasterYSize
    ulx = np.arange(0,size_x) * geo_out[1] + geo_out[0] + 0.5 *geo_out[1]
    uly = np.arange(0,size_y) * geo_out[5] + geo_out[3] + 0.5 *geo_out[5]   
    ulxx = np.zeros([len(uly),len(ulx)])
    ulyy = np.zeros([len(uly),len(ulx)])
    ulxx[:,:]=ulx[None,:]
    ulyy[:,:]=uly[:,None]
    
    ulx_flat = ulxx.flatten()
    uly_flat = ulyy.flatten()
    
    proj = osr.SpatialReference(wkt=dest.GetProjection())
    EPSG = int(proj.GetAttrValue('AUTHORITY',1))
    
    if EPSG != 4326:   
        inProj = Proj(init='epsg:%d' %EPSG)
        outProj = Proj(init='epsg:4326')
        ulx_flat, uly_flat = transform(inProj,outProj,ulx_flat, uly_flat)

    ulx_wgs = ulx_flat.reshape([len(uly),len(ulx)])              
    uly_wgs = uly_flat.reshape([len(uly),len(ulx)])            
           
    ulx_wgs_dis = np.abs(ulx_wgs - Coordinate[0])
    uly_wgs_dis = np.abs(uly_wgs - Coordinate[1])
    tot_dis = ulx_wgs_dis + uly_wgs_dis
    if np.nanmin(tot_dis) > (4 * (np.abs(ulx_wgs[0,0]-ulx_wgs[1,1]) + np.abs(uly_wgs[0,0]-uly_wgs[1,1]))):
        yID, xID = [np.nan, np.nan]
    else:
        yID, xID = np.argwhere(tot_dis == np.nanmin(tot_dis))[0]
    
    return(yID, xID)     
    
 
    
'''   
 
input_folders = [r"G:\SEBAL_Tadla\Daily_SEBAL\Dataset_out\Albedo_daily",r"G:\SEBAL_Tadla\Daily_SEBAL\Dataset_out\PreSEBAL_out\ALBEDO_MODIS_test\Albedo_new","G:\SEBAL_Tadla\Daily_SEBAL\Dataset_out\PreSEBAL_SEBAL_out\Albedo"]
input_formats = ["*{yyyy}_{doy}.tif", "*{yyyy}{mm:02d}{dd:02d}.tif", "*{yyyy}{mm:02d}{dd:02d}.tif"]
Coordinate = [-6.80689, 32.46206]
Startdate = "2016-10-01"
Enddate = "2017-09-30"

plt = Visualize_Graph(input_folders, input_formats, Coordinate, Startdate, Enddate)   
    
'''
    
    
    
    
    
    # Open data and get the value
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

