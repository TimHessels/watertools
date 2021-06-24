# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/DEM
"""

# General modules
import os
from ftplib import FTP
import urllib
import numpy as np
from osgeo import gdal
import time

import watertools.General.data_conversions as DC
import watertools.General.raster_conversions as RC

def DownloadData(output_folder, latlim, lonlim, dataset, level = None):
    """
    This function downloads SoilGrids data from SoilGrids.org

    Keyword arguments:
    output_folder -- directory of the result
	 latlim -- [ymin, ymax] (values must be between -50 and 50)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    level -- "sl1" .... "sl7"
    dataset -- ground dataset
    """
    
    dict_levels = dict()
    dict_levels["sl1"] = "0-5"
    dict_levels["sl2"] = "5-15"
    dict_levels["sl3"] = "15-30"
    dict_levels["sl4"] = "30-60"
    dict_levels["sl5"] = "60-100"
    dict_levels["sl6"] = "100-200"    
    
    if "conversion" in locals():
        del conversion
    
    # Define parameter depedent variables
    if dataset == "BULKDENSITY":
        nameEnd = os.path.join(output_folder, 'BulkDensity_%s_SoilGrids_kg-m-3.tif' %level)
        parameter = "bdod"
        conversion = 10 # cg/cm3 to kg/m3
        level_str = dict_levels[level]               
    if dataset == "NITROGEN":
        nameEnd = os.path.join(output_folder, 'Nitrogen_%s_SoilGrids_g_kg-1.tif' %level)        
        parameter = "nitrogen"   
        level_str = dict_levels[level]     
        conversion = 0.01         #cg/kg to g/kg        
    if dataset == "SOC":
        nameEnd = os.path.join(output_folder, 'SoilOrganicCarbonContent_%s_SoilGrids_g_kg.tif' %level)          
        parameter = "soc"
        level_str = dict_levels[level]      
        conversion = 0.1         #dg/kg to g/kg
    if dataset == "SOD":
        nameEnd = os.path.join(output_folder, 'SoilOrganicCarbonDensity_%s_SoilGrids_g_kg.tif' %level)        
        parameter = "ocd"
        conversion = 0.1          #dg/kg to g/kg
        level_str = dict_levels[level]             
    if dataset == "PH":
        nameEnd = os.path.join(output_folder, 'SoilPH_%s_SoilGrids_pH10.tif' %level) 
        parameter = "phh2o"       
        level_str = dict_levels[level]               
    if dataset == "CLAY":
        nameEnd = os.path.join(output_folder, 'ClayContentMassFraction_%s_SoilGrids_percentage.tif' %level)  
        parameter = "clay"   
        level_str = dict_levels[level]     
        conversion = 0.1          #g/kg to percentage         
    if dataset == "SAND":
        nameEnd = os.path.join(output_folder, 'SandContentMassFraction_%s_SoilGrids_percentage.tif' %level)        
        parameter = "sand"          
        level_str = dict_levels[level]    
        conversion = 0.1          #g/kg to percentage              
    if dataset == "SILT":
        nameEnd = os.path.join(output_folder, 'SiltContentMassFraction_%s_SoilGrids_percentage.tif' %level)     
        parameter = "silt"     
        level_str = dict_levels[level]    
        conversion = 0.1          #g/kg to percentage              

    if not os.path.exists(nameEnd):
            
        # Download, extract, and converts all the files to tiff files
        try:
                    
            url = "https://maps.isric.org/mapserv?map=/map/%s.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=%s_%scm_mean&FORMAT=image/tiff&SUBSET=long(%f,%f)&SUBSET=lat(%f,%f)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"  %(parameter, parameter, level_str, lonlim[0], lonlim[1], latlim[0], latlim[1])      
            #url = "http://85.214.241.121:8080/geoserver/ows?service=WCS&version=2.0.1&request=GetCoverage&CoverageId=%s_M_%s250m&subset=Long(%f,%f)&subset=Lat(%f,%f)" %(dataset, level_name, lonlim[0], lonlim[1], latlim[0], latlim[1])
            print(url)
            urllib.request.urlretrieve(url, filename=nameEnd)
            
            if "conversion" in locals():
                print("Conversion is applied of %s" %conversion)
                dest = gdal.Open(nameEnd)
                geo = dest.GetGeoTransform()
                proj = "WGS84"
                Array = dest.GetRasterBand(1).ReadAsArray()
                del dest
                time.sleep(1)
                Array = np.float_(Array) * conversion
                
                try:
                    Array = RC.gap_filling(Array, 0)
                except:
                    pass
                
                DC.Save_as_tiff(nameEnd, Array, geo, proj)

        except:
            print("Was not able to create the wanted dataset")
        
    return()


    

