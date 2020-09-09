# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 10:10:43 2020

@author: timhe
"""
import sys
from watertools.Collect.WAPOR.DataAccess import WAPOR

def main(output_folder, Startdate, Enddate, latlim, lonlim, Parameter, Area = None, Version = "2"):
    
    """
    This function downloads WAPOR data for the specified time
    interval, and spatial extent.

    Keyword arguments:
    -- output_folder: (string) to define the folder where to save the created tiff files.
    -- Startdate: (string) defines the startdate of the required dataset in the following format "yyyy-mm-dd"
    -- Enddate: (string) defines the enddate of the required dataset in the following format "yyyy-mm-dd"
    -- Latlim: (array) defines the latitude limits in the following format [Latitude_minimum, Latitude_maximum] e.g. [10, 13]
    -- Lonlim: (array) defines the longitude limits in the following format [Longitude_minimum, Longitude_maximum] e.g. [10, 13]
    -- API_KEY: (string) this is a private API KEY that can be collected here: https://io.apps.fao.org/gismgr/api/v1/swagger-ui.html#/IAM/apiKeySignIn
    -- Parameter: (string) define the parameter that is required (possible Parameters: watertools.Collect.WAPOR.DataAccess.VariablesInfo().descriptions.keys())
    -- version: (string) default is 2, but if version 1 is required use this parameter (options are "1" or "2")

    """
    
    
    print('\nDownload WAPOR %s data for period %s till %s' %(Parameter, Startdate, Enddate))
    WAPOR(output_folder, Startdate, Enddate, latlim, lonlim, Parameter, Area = None, Version = "2")

if __name__ == '__main__':
    main(sys.argv)   