# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 14:11:15 2018

@author: tih
"""
import numpy as np
import pandas as pd
import glob
import os
import gdal

import watertools.General.raster_conversions as RC

def Visualize_Graph(input_folders, input_formats, Startdate, Enddate, input_shp, obj = 'id', Shp_value = 1, vmin = None, vmax = None):

    Dates = pd.date_range(Startdate, Enddate, freq = "D")
    
    Datasets_end = np.ones([len(Dates), len(input_folders)+1]) * np.nan
    
    for j in range(0,len(input_folders)):
    
        input_folder = input_folders[j]  
        input_format = input_formats[j]
        Dataset = Get_Dataset_Polygon(input_folder, input_format, Dates, input_shp, Shp_prop = obj, Shp_value = Shp_value)   
        
        if j == 0:
            
            Datasets_end[:,0:2] = Dataset
            
        else:
            Datasets_end[:,j + 1] = Dataset[:,1]
    
    import matplotlib.pyplot as plt
    
    colors = ['red','blue','green','black','yellow']
    
    for k in range(0,len(input_folders)):
        plt.scatter(Dates, Datasets_end[:,k+1], s =10, c = colors[k])
    plt.ylim(vmin, vmax)
    plt.legend(input_folders, loc='upper center',bbox_to_anchor=(0.5, -0.1))
    plt.show()  
    
    return(plt)
       
def Get_Dataset_Polygon(input_folder, input_format, Dates, input_shp, Shp_prop = 'id', Shp_value = 1, Stats = "mean"):
      
    Dataset = np.ones([len(Dates), 2]) * np.nan
    i = 0
    os.chdir(input_folder)
    
    filename_test = glob.glob('*.tif')[0]
    filename_in = os.path.join(input_folder, filename_test)
    dest = gdal.Open(filename_in) 
    geo = dest.GetGeoTransform()  
    
    epsg_tiff = RC.Get_epsg(dest)
    epsg_shp = RC.Get_epsg(input_shp, 'shp')
    
    if epsg_tiff != epsg_shp:
        input_shp = RC.reproject_shapefile(input_shp, epsg_tiff)
        
    # Create Mask of polygon
    mask_tiff = RC.GDAL_rasterize(input_shp, geo[1], Shp_prop)
    
    dest = RC.reproject_dataset_example(mask_tiff, filename_in)
           
    MASK = dest.GetRasterBand(1).ReadAsArray()
    MASK[MASK!=Shp_value] = np.nan           
    MASK[MASK==Shp_value] = 1

    for Date in Dates:

        Day = Date.day
        Month = Date.month
        Year = Date.year
        DOY = Date.dayofyear
        
        filename = glob.glob(input_format.format(yyyy=Year, mm=Month, dd=Day, doy = DOY))
        
        Dataset[i,0] = Date.toordinal()
        
    
        if len(filename) > 0: 
            filename = filename[0]
           
            
            dest = gdal.Open(filename)
            Array = dest.GetRasterBand(1).ReadAsArray()
            Array = np.float_(Array)
            Array[Array==-9999] = np.nan
            try:
                if Stats == "mean":
                    Value = np.nanmean(Array[MASK==1])
                if Stats == "std":
                    Value = np.nanstd(Array[MASK==1])   
                
            except:
                dest = RC.reproject_dataset_example(mask_tiff, filename)
           
                MASK = dest.GetRasterBand(1).ReadAsArray()
                MASK[MASK>0] = 1
                MASK[MASK<=0] = np.nan   
                if Stats == "mean":
                    Value = np.nanmean(Array[MASK==1])
                if Stats == "std":
                    Value = np.nanstd(Array[MASK==1])                
                
                
            Dataset[i,1] = Value
        i += 1
        
    return(Dataset)