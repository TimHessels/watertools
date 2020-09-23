"""
Description:
This tool collects FAO's WAPOR data of LEVEL 1 and LEVEL 2 (excluding the seasonal paramters). The output will be tiff files and no conversions is 
preformed to the units. The output are in the same units as on the WAPOR Portal (https://wapor.apps.fao.org/home/1)

The following input is required:
    
    -- output_folder: (string) to define the folder where to save the created tiff files.
    -- Startdate: (string) defines the startdate of the required dataset in the following format "yyyy-mm-dd"
    -- Enddate: (string) defines the enddate of the required dataset in the following format "yyyy-mm-dd"
    -- Latlim: (array) defines the latitude limits in the following format [Latitude_minimum, Latitude_maximum] e.g. [10, 13]
    -- Lonlim: (array) defines the longitude limits in the following format [Longitude_minimum, Longitude_maximum] e.g. [10, 13]
    -- API_KEY: (string) this is a private API KEY that can be collected here: https://io.apps.fao.org/gismgr/api/v1/swagger-ui.html#/IAM/apiKeySignIn
    -- Parameter: (string) define the parameter that is required
    -- version: (string) default is 2, but if version 1 is required use this parameter (options are "1" or "2")
    
To get a describtion of all the possible parameters type:
    
import WaporAPI
WaporAPI.Collect.VariablesInfo.descriptions 

Examples:
import watertools
watertools.Collect.WAPOR.Get_Layer(r"G:\Project_MetaMeta\Input_Data2", "2011-01-01", "2011-02-01", [8.2, 8.7], [39, 39.5], "L2_T_D")
"""

from .Get_Layer import main as Get_Layer


__all__ = ['Get_Layer']
__version__ = '0.0.0'