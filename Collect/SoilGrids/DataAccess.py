# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/DEM
"""

# General modules
import os
from ftplib import FTP
import urllib

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
    # Define parameter depedent variables
    if dataset == "BLDFIE":
        nameEnd = os.path.join(output_folder, 'BulkDensity_%s_SoilGrids_kg-m-3.tif' %level)
        level_name = "%s_" %level
    if dataset == "CLYPPT":
        nameEnd = os.path.join(output_folder, 'ClayContentMassFraction_%s_SoilGrids_percentage.tif' %level)
        level_name = "%s_" %level        
    if dataset == "ORCDRC":
        nameEnd = os.path.join(output_folder, 'SoilOrganicCarbonContent_%s_SoilGrids_g_kg.tif' %level)        
        level_name = "%s_" %level
    if dataset == "OCSTHA":
        nameEnd = os.path.join(output_folder, 'SoilOrganicCarbonStock_%s_SoilGrids_tonnes-ha-1.tif' %level)        
        level_name = "%s_" %level
    if dataset == "CRFVOL":
        nameEnd = os.path.join(output_folder, 'CoarseFragmentVolumetric_%s_SoilGrids_percentage.tif' %level)        
        level_name = "%s_" %level
    if dataset == "SLTPPT":
        nameEnd = os.path.join(output_folder, 'SiltContentMassFraction_%s_SoilGrids_percentage.tif' %level)        
        level_name = "%s_" %level
    if dataset == "SNDPPT":
        nameEnd = os.path.join(output_folder, 'SandContentMassFraction_%s_SoilGrids_percentage.tif' %level)     
        level_name = "%s_" %level
    if dataset == "BDRICM":
        nameEnd = os.path.join(output_folder, 'DepthToBedrock_SoilGrids_cm.tif') 
        level_name = ""
    if dataset == "BDTICM":
        nameEnd = os.path.join(output_folder, 'AbsoluteDepthToBedrock_SoilGrids_cm.tif') 
        level_name = ""
    if dataset == "BDRLOG":
        nameEnd = os.path.join(output_folder, 'PredictedProbabilityOfOccurrence_SoilGrids_percentage.tif') 
        level_name = ""
    if dataset == "PHIKCL":
        nameEnd = os.path.join(output_folder, 'SoilPH_%s_SoilGrids_KCi10.tif' %level) 
        level_name = "%s_" %level


    if not os.path.exists(nameEnd):
            
        # Download, extract, and converts all the files to tiff files
        try:
                    
            url = "http://85.214.241.121:8080/geoserver/ows?service=WCS&version=2.0.1&request=GetCoverage&CoverageId=%s_M_%s250m&subset=Long(%f,%f)&subset=Lat(%f,%f)" %(dataset, level_name, lonlim[0], lonlim[1], latlim[0], latlim[1])
            #print(url)
            urllib.request.urlretrieve(url, filename=nameEnd)

        except:
            print("Was not able to create the wanted dataset")
        
    return()


    

