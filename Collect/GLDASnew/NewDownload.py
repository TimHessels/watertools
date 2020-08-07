# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 08:19:31 2020

@author: timhe
"""
import os
import requests
import datetime
import pandas as pd

parameters = ["Psurf_f_inst"]
latlim = [25.313, 59.063]
lonlim = [-121.641, -82.969]
Startdate = "2001-01-01"
Enddate = "2001-01-03"
Version = "2.1"
Model = "GLDAS_NOAH025_M"
output_folder = r"F:\PHD\PHD_2\Datasets_Mead_Farm"

"GLDAS_NOAH025_M"  2.1 en 2.0
"GLDAS_NOAH025_3H" 2.1 en 2.0
"GLDAS_NOAH10_3H"  2.1 en 2.0
"GLDAS_NOAH10_M"   2.1 en 2.0

"GLDAS_CLSM025_D"  2.0
"GLDAS_CLSM10_3H"  2.1
"GLDAS_CLSM10_M"   2.1

"GLDAS_VIC10_M"    2.1
"GLDAS_VIC10_3H"   2.1








dates = pd.date_range(Startdate, Enddate, freq = "3H")
filename_temp_nc_format = "temp_{model}_{version}_{year}{month:02d}{day:02d}_{hour:02d}{minute:02d}_{parameter_list}.nc"

for date in dates:
    
    year = date.year
    month = date.month
    day = date.day
    doy = date.dayofyear
    hour = date.hour
    minute = date.minute
    version = Version.replace(".","")
    version_dot = Version
    parameter_list = "-".join(parameters)
    parameter = "%2C".join(parameters)
    url_format = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2F{model_first}%2F{model}.{version_dot}%2F{year}%2F{model}.A{year}{month:02d}.0{version}.nc4&FORMAT=bmM0Lw&BBOX={latmin}%2C{lonmin}%2C{latmax}%2C{lonmax}&LABEL={model}.A{year}{month:02d}.0{version}.nc4.SUB.nc4&SHORTNAME={model}&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION={version_dot}&VARIABLES={para}"

    url = url_format.format(year = year, month = month, latmin = latlim[0], latmax = latlim[1], lonmin = lonlim[0], lonmax = lonlim[1], para = parameter, version = version, version_dot = version_dot, model = Model, model_first = Model.split("_")[0])
    
    filename_temp_nc = filename_temp_nc_format.format(year = year, month = month, day = day, model = Model, version = version, hour = hour, minute = minute, parameter_list = parameter_list)
    download_file(url, os.path.join(output_folder, filename_temp_nc))




def download_file(dl_url, local_filename):
    
    session = requests.Session()
    downloaded = 0
    with session.get(dl_url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    downloaded += len(chunk)
                    # f.flush()
    return 


class VariablesInfo_NOAH:
   parameters = {'albedo_inst': 'Albedo_inst',
                'avgsurft_inst': 'AvgSurfT_inst',
                'canopint_inst': 'CanopInt_inst',
                'ecanop_tavg': 'ECanop_tavg',
                'esoil_tavg': 'ESoil_tavg',
                'evap_tavg': 'Evap_tavg',
                'lwdown_f_tavg': 'LWdown_f_tavg',
                'lwnet_tavg': 'Lwnet_tavg',
                'potevap_tavg': 'PotEvap_tavg',
                'psurf_f_inst': 'Psurf_f_inst',
                'qair_f_inst': 'Qair_f_inst',
                'qg_tavg': 'Qg_tavg',
                'qh_tavg': 'Qh_tavg',
                'qle_tavg': 'Qle_tavg',
                'qs_acc': 'Qs_acc',
                'qsb_acc': 'Qsb_acc',
                'qsm_acc': 'Qsm_acc',
                'rainf_f_tavg': 'Rainf_f_tavg',
                'rainf_tavg': 'Rainf_tavg',
                'rootmoist_inst': 'RootMoist_inst',
                'snowdepth_inst': 'SnowDepth_inst',
                'snowf_tavg': 'Snowf_tavg',
                'sm0_10cm_ins': 'SoilMoi0_10cm_inst',
                'sm10_40cm_ins': 'SoilMoi10_40cm_inst',
                'sm40_100cm_ins': 'SoilMoi40_100cm_inst',
                'sm100_200cm_ins': 'SoilMoi100_200cm_inst',
                'st0_10cm_ins': 'SoilTMP0_10cm_inst',
                'st10_40cm_ins': 'SoilTMP10_40cm_inst',
                'st40_100cm_ins': 'SoilTMP40_100cm_inst',
                'st100_200cm_ins':'SoilTMP100_200cm_inst',
                'swdown_f_tavg': 'SWdown_f_tavg',
                'swe_inst': 'SWE_inst',
                'swnet_tavg': 'Swnet_tavg',
                'tair_f_inst': 'Tair_f_inst',
                'tveg_tavg' : 'Tveg_tavg',
                'wind_f_inst': 'Wind_f_inst'}
   names = {'albedo_inst': 'Albedo',
             'avgsurft_inst': 'T',
             'canopint_inst': 'TotCanopyWaterStorage',
             'ecanop_tavg': 'ECanop',
             'esoil_tavg': 'ESoil',             
             'evap_tavg': 'ET',
             'lwdown_f_tavg': 'LWdown',
             'lwnet_tavg': 'LWnet',
             'potevap_tavg': 'Epot'             
             'psurf_f_inst': 'P',
             'qair_f_inst': 'Hum',
             'qg_tavg': 'G',
             'qh_tavg': 'H',
             'qle_tavg': 'LE',
             'qs_acc': 'Rsur',
             'qsb_acc': 'Rsubsur',
             'qsm_acc': 'SnowMelt',
             'rainf_f_tavg': 'P',
             'rainf_tavg': 'Prain',
             'rootmoist_inst': 'RootMoist',
             'snowdepth_inst': 'SnowDepth',           
             'swe_inst': 'SnowWaterEquivalent',
             'swdown_f_tavg': 'SWdown',
             'swnet_tavg': 'SWnet',
             'snowf_tavg': 'Snow',
             'sm0_10cm_ins': 'S1',
             'sm10_40cm_ins': 'S2',
             'sm40_100cm_ins': 'S3',
             'sm100_200cm_ins': 'S4',
             'st0_10cm_ins': 'Ts1',
             'st10_40cm_ins': 'Ts2',
             'st40_100cm_ins': 'Ts3',
             'st100_200cm_ins': 'Ts4',
             'tair_f_inst': 'Tair',
             'wind_f_inst': 'W',
             'tveg_tavg' : 'Transpiration'}
    descriptions = {'albedo_inst': 'Albedo [percentage]',
                    'avgsurft_inst': 'surface average skin surface temperature [k]',
                    'canopint_inst': 'surface plant canopy surface water [kg/m^2]',
                    'ecanop_tavg': 'Canopy water evaporation [W/m^2]',
                    'esoil_tavg': 'Direct evaporation from bare soil [W/m^2]',    
                    'evap_tavg': 'surface total evapotranspiration [kg/m^2/s]',
                    'lwdown_f_tavg': 'surface surface incident longwave radiation'
                               ' [w/m^2]',
                    'lwnet_tavg': 'surface net longwave radiation [w/m^2]',
                    'potevap_tavg': 'Potential evaporation rate [W/m^2]',
                    'psurf_f_inst': 'surface surface pressure [kPa]',
                    'qair_f_inst': 'surface near surface specific humidity [kg/kg]',
                    'qg_tavg': 'surface ground heat flux [w/m^2]',
                    'qh_tavg': 'surface sensible heat flux [w/m^2]',
                    'qle_tavg': 'surface latent heat flux [w/m^2]',
                    'qs_acc': 'storm surface runoff [kg/m^2/s]',
                    'qsb_acc': 'baseflow-groundwater runoff [kg/m^2/s]',
                    'qsm_acc': 'surface snowmelt [kg/m^2/s]',
                    'rainf_f_tavg': 'Total precipitation rate [kg/m^2/s]',
                    'rainf_tavg': 'Rain precipitation rate [kg/m^2/s]',
                    'rootmoist_inst': 'Root zone soil moisture [kg/m^2]',
                    'snowdepth_inst': 'Snow depth [m]',   
                    'swe_inst': 'surface snow water equivalent [kg/m^2]',
                    'swdown_f_tavg': 'surface surface incident shortwave radiation'
                               ' [w/m^2]',
                    'swnet_tavg': 'surface net shortwave radiation [w/m^2]',
                    'snowf_tavg': 'surface snowfall rate [kg/m^2/s]',
                    'sm0_10cm_ins': '0-10 cm underground soil moisture content'
                               ' [kg/m^2]',
                    'sm10_40cm_ins': '10-40 cm underground soil moisture content'
                               ' [kg/m^2]',
                    'sm40_100cm_ins': '40-100 cm underground soil moisture content'
                               ' [kg/m^2]',
                    'sm100_200cm_ins': '100-200 cm underground soil moisture content'
                               ' [kg/m^2]',
                    'st0_10cm_ins': '0-10 cm underground soil temperature [k]',
                    'st10_40cm_ins': '10-40 cm underground soil temperature [k]',
                    'st40_100cm_ins': '40-100 cm underground soil temperature [k]',
                    'st100_200cm_ins': '100-200 cm underground soil temperature [k]',
                    'tair_f_inst': 'surface near surface air temperature [k]',
                    'wind_f_inst': 'surface near surface wind speed [m/s]',
                    'tveg_tavg' : 'transpiration [w/m^2]'}
    factors = {'albedo_inst': 1,
               'avgsurft_inst': -273.15,
               'canopint_inst': 1,
               'ecanop_tavg': 1,
               'esoil_tavg': 1,                  
               'evap_tavg': 86400,# !!! delen door uren/24
               'lwdown_f_tavg': 1,
               'lwnet_tavg': 1,
               'potevap_tavg': 1,
               'psurf_f_inst': 0.001,
               'qair_f_inst': 1,
               'qg_tavg': 1,
               'qh_tavg': 1,
               'qle_tavg': 1,
               'qs_acc': 86400,# !!! delen door uren/24
               'qsb_acc': 86400,# !!! delen door uren/24
               'qsm_acc': 86400,# !!! delen door uren/24
               'rainf_f_tavg': 86400, # !!! delen door uren/24
               'rainf_tavg': 86400,               
               'rootmoist_inst': 1,
               'snowdepth_inst': 1,   
               'swe_inst': 1,
               'swdown_f_tavg': 1,
               'swnet_tavg': 1,
               'snowf_tavg': 1,
               'sm0_10cm_ins': 1,
               'sm10_40cm_ins': 1,
               'sm40_100cm_ins': 1,
               'sm100_200cm_ins': 1,
               'st0_10cm_ins': -273.15,
               'st10_40cm_ins': -273.15,
               'st40_100cm_ins': -273.15,
               'st100_200cm_ins': -273.15,
               'tair_f_inst': -273.15,
               'wind_f_inst': 1,
               'tveg_tavg' : 1}
    types = {'albedo_inst': 'state',
             'avgsurft_inst': 'state',
             'canopint_inst': 'state',
             'ecanop_tavg': 'state',
             'esoil_tavg': 'state',   
             'evap_tavg': 'flux',
             'lwdown_f_tavg': 'state',
             'lwnet_tavg': 'state',
             'potevap_tavg': 'state',             
             'psurf_f_inst': 'state',
             'qair_f_inst': 'state',
             'qg_tavg': 'state',
             'qh_tavg': 'state',
             'qle_tavg': 'state',
             'qs_acc': 'flux',
             'qsb_acc': 'flux',
             'qsm_acc': 'flux',
             'rainf_f_tavg': 'flux',
             'rainf_tavg': 'flux',    
             'rootmoist_inst': 'state',
             'snowdepth_inst': 'state',                
             'swe_inst': 'state',
             'swdown_f_tavg': 'state',
             'swnet_tavg': 'state',
             'snowf_tavg': 'state',
             'sm0_10cm_ins': 'state',
             'sm10_40cm_ins': 'state',
             'sm40_100cm_ins': 'state',
             'sm100_200cm_ins': 'state',
             'st0_10cm_ins': 'state',
             'st10_40cm_ins': 'state',
             'st40_100cm_ins': 'state',
             'st100_200cm_ins': 'state',
             'tair_f_inst': 'state',
             'wind_f_inst': 'state',
             'tveg_tavg' : 'state'}

    def __init__(self, step):
        if step == 'three_hourly':
            self.units = {'albedo_inst': 'Percentage',
                          'avgsurft_inst': 'C',
                          'canopint_inst': 'mm',
                          'ecanop_tavg': 'W-m-2',
                          'esoil_tavg': 'W-m-2',                          
                          'evap_tavg': 'mm-3hour-1',
                          'lwdown_f_tavg': 'W-m-2',
                          'lwnet_tavg': 'W-m-2',
                          'potevap_tavg': 'W-m-2',   
                          'psurf_f_inst': 'kpa',
                          'qair_f_inst': 'kg-kg',
                          'qg_tavg': 'W-m-2',
                          'qh_tavg': 'W-m-2',
                          'qle_tavg': 'W-m-2',
                          'qs_acc': 'mm-3hour-1',
                          'qsb_acc': 'mm-3hour-1',
                          'qsm_acc': 'mm-3hour-1',
                          'rainf_f_tavg': 'mm-3hour-1',
                          'rainf_tavg':  'mm-3hour-1',   
                          'rootmoist_inst': 'mm',
                          'snowdepth_inst': 'm',                             
                          'swe_inst': 'mm',
                          'swdown_f_tavg': 'W-m-2',
                          'swnet_tavg': 'W-m-2',
                          'snowf_tavg': 'mm',
                          'sm0_10cm_ins': 'kg-m-2',
                          'sm10_40cm_ins': 'kg-m-2',
                          'sm40_100cm_ins': 'kg-m-2',
                          'sm100_200cm_ins': 'kg-m-2',
                          'st0_10cm_ins': 'C',
                          'st10_40cm_ins': 'C',
                          'st40_100cm_ins': 'C',
                          'st100_200cm_ins': 'C',
                          'tair_f_inst': 'C',
                          'wind_f_inst': 'm-s-1',
                          'tveg_tavg' : 'W-m-2'}
        elif step == 'daily':
            self.units = {'albedo_inst': 'Percentage',
                          'avgsurft_inst': 'C',
                          'canopint_inst': 'mm',
                          'ecanop_tavg': 'W-m-2',
                          'esoil_tavg': 'W-m-2',                           
                          'evap_tavg': 'mm-day-1',
                          'lwdown_f_tavg': 'W-m-2',
                          'lwnet_tavg': 'W-m-2',
                          'potevap_tavg': 'W-m-2',   
                          'psurf_f_inst': 'kpa',
                          'qair_f_inst': 'kg-kg',
                          'qg_tavg': 'W-m-2',
                          'qh_tavg': 'W-m-2',
                          'qle_tavg': 'W-m-2',
                          'qs_acc': 'mm-day-1',
                          'qsb_acc': 'mm-day-1',
                          'qsm_acc': 'mm-day-1',
                          'rainf_f_tavg': 'mm-day-1',
                          'rainf_tavg':  'mm-day-1',                             
                          'rootmoist_inst': 'mm',
                          'snowdepth_inst': 'm',  
                          'swe_inst': 'mm',
                          'swdown_f_tavg': 'W-m-2',
                          'swnet_tavg': 'W-m-2',
                          'snowf_tavg': 'mm',
                          'sm0_10cm_ins': 'kg-m-2',
                          'sm10_40cm_ins': 'kg-m-2',
                          'sm40_100cm_ins': 'kg-m-2',
                          'sm100_200cm_ins': 'kg-m-2',
                          'st0_10cm_ins': 'C',
                          'st10_40cm_ins': 'C',
                          'st40_100cm_ins': 'C',
                          'st100_200cm_ins': 'C',
                          'tair_f_inst': 'C',
                          'wind_f_inst': 'm-s-1',
                          'tveg_tavg' : 'W-m-2'}
        elif step == 'monthly':
            self.units = {'albedo_inst': 'Percentage',
                          'avgsurft_inst': 'C',
                          'canopint_inst': 'mm',
                          'ecanop_tavg': 'W-m-2',
                          'esoil_tavg': 'W-m-2',                           
                          'evap_tavg': 'mm-month-1',
                          'lwdown_f_tavg': 'W-m-2',
                          'lwnet_tavg': 'W-m-2',
                          'potevap_tavg': 'W-m-2', 
                          'psurf_f_inst': 'kpa',
                          'qair_f_inst': 'kg-kg',
                          'qg_tavg': 'W-m-2',
                          'qh_tavg': 'W-m-2',
                          'qle_tavg': 'W-m-2',
                          'qs_acc': 'mm-month-1',
                          'qsb_acc': 'mm-month-1',
                          'qsm_acc': 'mm-month-1',
                          'rainf_f_tavg': 'mm-month-1',
                          'rainf_tavg':  'mm-month-1', 
                          'rootmoist_inst': 'mm',
                          'snowdepth_inst': 'm',                            
                          'swe_inst': 'mm',
                          'swdown_f_tavg': 'W-m-2',
                          'swnet_tavg': 'W-m-2',
                          'snowf_tavg': 'mm',
                          'sm0_10cm_ins': 'kg-m-2',
                          'sm10_40cm_ins': 'kg-m-2',
                          'sm40_100cm_ins': 'kg-m-2',
                          'sm100_200cm_ins': 'kg-m-2',
                          'st0_10cm_ins': 'C',
                          'st10_40cm_ins': 'C',
                          'st40_100cm_ins': 'C',
                          'st100_200cm_ins': 'C',
                          'tair_f_inst': 'C',
                          'wind_f_inst': 'm-s-1',
                          'tveg_tavg' : 'W-m-2'}
        else:
            raise KeyError("The input time step is not supported")


class VariablesInfo_CLSM:
    
    """
    This class contains the information about the GLDAS variables
    """
ACond_tavg = Aerodynamic conductance (m s-1)
 = Average surface skin temperature (K)
 = Plant canopy surface water (kg m-2)
ECanop_tavg = Canopy water evaporation (kg m-2 s-1)
ESoil_tavg = Direct evaporation from bare soil (kg m-2 s-1)
 = Evapotranspiration (kg m-2 s-1)
EvapSnow_tavg = Snow evaporation (kg m-2 s-1)
GWS_tavg = Ground water storage (mm)
 = Downward longwave radiation flux (W m-2)
 = Net longwave radiation flux (W m-2)
 = Surface air pressure (Pa)
 = Specific humidity (kg kg-1)
 = Ground heat flux (W m-2)
 = Sensible heat net flux (W m-2)
 = Latent heat net flux (W m-2)
 = Storm surface runoff (kg m-2 s-1)
 = Baseflow-groundwater runoff (kg m-2 s-1)
 = Snow melt (kg m-2 s-1)
 = Total precipitation rate (kg m-2 s-1)
 = Rain precipitation rate (kg m-2 s-1)
SnowDepth_tavg = Snow depth (m)
 = Snow precipitation rate (kg m-2 s-1)
SnowT_tavg = Snow surface temperature (K)
 = Profile soil moisture (kg m-2)
 = Root zone soil moisture (kg m-2)
 = Surface soil moisture (kg m-2)
 = Downward shortwave radiation flux (W m-2)
 = Snow depth water equivalent (kg m-2)
 = Net shortwave radiation flux (W m-2)
 = Air temperature (K)
 = Transpiration (kg m-2 s-1)
TWS_tavg = Terrestrial water storage (mm)
 = Wind speed (m s-1)

    parameters = {'avgsurft_tavg': 'AvgSurfT_tavg',
             'canopint_tavg': 'CanopInt_tavg',
             'evap_tavg': 'Evap_tavg',
             'lwdown_f_tavg': 'LWdown_f_tavg',
             'lwnet_tavg': 'Lwnet_tavg',
             'psurf_f_tavg': 'Psurf_f_tavg',
             'qair_f_tavg': 'Qair_f_tavg',
             'qg_tavg': 'Qg_tavg',
             'qh_tavg': 'Qh_tavg',
             'qle_tavg': 'Qle_tavg',
             'qs_tavg': 'Qs_tavg',
             'qsb_tavg': 'Qsb_tavg',
             'qsm_tavg': 'Qsm_tavg',
             'rainf_f_tavg': 'Rainf_f_tavg',
             'rainf_tavg': 'Rainf_tavg',
             'swe_tavg': 'SWE_tavg',
             'swdown_f_tavg': 'SWdown_f_tavg',
             'swnet_tavg': 'SWnet_tavg',
             'snowf_tavg': 'Snowf_tavg',
             'soilmoist_s_tav': 'SoilMoist_S_tavg',
             'soilmoist_rz_ta': 'SoilMoist_RZ_tavg',
             'soilmoist_p_tav': 'SoilMoist_P_tavg',
             'tair_f_tavg': 'Tair_f_tavg',
             'wind_f_tavg': 'Wind_f_tavg',
             'tveg_tavg' : 'TVeg_tavg'}  
    
    
    names = {'avgsurft_tavg': 'SurfaceTemperature',
             'canopint_tavg': 'TotCanopyWaterStorage',
             'evap_tavg': 'ET',
             'lwdown_f_tavg': 'LWdown',
             'lwnet_tavg': 'LWnet',
             'psurf_f_tavg': 'P',
             'qair_f_tavg': 'Hum',
             'qg_tavg': 'G',
             'qh_tavg': 'H',
             'qle_tavg': 'LE',
             'qs_tavg': 'Rsur',
             'qsb_tavg': 'Rsubsur',
             'qsm_tavg': 'SnowMelt',
             'rainf_f_tavg': 'P',
             'rainf_tavg': 'Prain',
             'swe_tavg': 'SnowWaterEquivalent',
             'swdown_f_tavg': 'SWdown',
             'swnet_tavg': 'SWnet',
             'snowf_tavg': 'Snow',
             'soilmoist_s_tav': 'SoilMoisturSurface',
             'soilmoist_rz_ta': 'SoilMoistureRootZone',
             'soilmoist_p_tav': 'SoilMoistureProfile',
             'tair_f_tavg': 'Tair',
             'wind_f_tavg': 'W',
             'tveg_tavg' : 'Transpiration'}
    
    descriptions = {'avgsurft_tavg': 'surface average surface temperature [k]',
                    'canopint_tavg': 'surface plant canopy surface water [kg/m^2]',
                    'evap_tavg': 'surface total evapotranspiration [kg/m^2/s]',
                    'lwdown_f_tavg': 'surface surface incident longwave radiation'
                               ' [w/m^2]',
                    'lwnet_tavg': 'surface net longwave radiation [w/m^2]',
                    'psurf_f_tavg': 'surface surface pressure [kPa]',
                    'qair_f_tavg': 'surface near surface specific humidity [kg/kg]',
                    'qg_tavg': 'surface ground heat flux [w/m^2]',
                    'qh_tavg': 'surface sensible heat flux [w/m^2]',
                    'qle_tavg': 'surface latent heat flux [w/m^2]',
                    'qs_tavg': 'storm surface runoff [kg/m^2/s]',
                    'qsb_tavg': 'baseflow-groundwater runoff [kg/m^2/s]',
                    'qsm_tavg': 'surface snowmelt [kg/m^2/s]',
                    'rainf_f_tavg': 'surface rainfall rate [kg/m^2/s]',
                    'rainf_tavg': 'Rain precipitation rate [kg/m^2/s]',
                    'swe_tavg': 'surface snow water equivalent [kg/m^2]',
                    'swdown_f_tavg': 'surface surface incident shortwave radiation'
                               ' [w/m^2]',
                    'swnet_tavg': 'surface net shortwave radiation [w/m^2]',
                    'snowf_tavg': 'surface snowfall rate [kg/m^2/s]',
                    'soilmoist_s_tav': 'surface soil moisture [kg/m^2]',
                    'soilmoist_rz_ta': 'root zone soil moisture [kg/m^2]',
                    'soilmoist_p_tav': 'profile soil moisture [kg/m^2]',
                    'tair_f_tavg': 'surface near surface air temperature [k]',
                    'wind_f_tavg': 'surface near surface wind speed [m/s]',
                    'tveg_tavg' : 'transpiration [w/m^2]'}
    
    factors = {'avgsurft_tavg': -273.15,
               'canopint_tavg': 1,
               'evap_tavg': 86400,
               'lwdown_f_tavg': 1,
               'lwnet_tavg': 1,
               'psurf_f_tavg': 0.001,
               'qair_f_tavg': 1,
               'qg_tavg': 1,
               'qh_tavg': 1,
               'qle_tavg': 1,
               'qs_tavg': 86400,
               'qsb_tavg': 86400,
               'qsm_tavg': 86400,
               'rainf_f_tavg': 86400,
               'rainf_tavg': 86400,
               'swe_tavg': 1,
               'swdown_f_tavg': 1,
               'swnet_tavg': 1,
               'snowf_tavg': 1,
               'soilmoist_s_tav': 1,
               'soilmoist_rz_ta': 1,
               'soilmoist_p_tav': 1,
               'tair_f_tavg': -273.15,
               'wind_f_tavg': 0.75,
               'tveg_tavg' : 1}
    types = {'avgsurft_tavg': 'state',
             'canopint_tavg': 'state',
             'evap_tavg': 'flux',
             'lwdown_f_tavg': 'state',
             'lwnet_tavg': 'state',
             'psurf_f_tavg': 'state',
             'qair_f_tavg': 'state',
             'qg_tavg': 'state',
             'qh_tavg': 'state',
             'qle_tavg': 'state',
             'qs_tavg': 'flux',
             'qsb_tavg': 'flux',
             'qsm_tavg': 'flux',
             'rainf_f_tavg': 'flux',
             'rainf_tavg': 'flux',
             'swe_tavg': 'state',
             'swdown_f_tavg': 'state',
             'swnet_tavg': 'state',
             'snowf_tavg': 'state',
             'soilmoist_s_tav': 'state',
             'soilmoist_rz_ta': 'state',
             'soilmoist_p_tav': 'state',
             'tair_f_tavg': 'state',
             'wind_f_tavg': 'state',
             'tveg_tavg' : 'state'}

    def __init__(self, step):
        if step == 'three_hourly':
            self.units = {'avgsurft_tavg': 'C',
                          'canopint_tavg': 'mm',
                          'evap_tavg': 'mm-3hour-1',
                          'lwdown_f_tavg': 'W-m-2',
                          'lwnet_tavg': 'W-m-2',
                          'psurf_f_tavg': 'kpa',
                          'qair_f_tavg': 'kg-kg',
                          'qg_tavg': 'W-m-2',
                          'qh_tavg': 'W-m-2',
                          'qle_tavg': 'W-m-2',
                          'qs_tavg': 'mm-3hour-1',
                          'qsb_tavg': 'mm-3hour-1',
                          'qsm_tavg': 'mm-3hour-1',
                          'rainf_f_tavg': 'mm-3hour-1',
                          'rainf_tavg': 'mm-3hour-1',
                          'swe_tavg': 'mm',
                          'swdown_f_tavg': 'W-m-2',
                          'swnet_tavg': 'W-m-2',
                          'snowf_tavg': 'mm',
                          'soilmoist_s_tav': 'mm',
                          'soilmoist_rz_ta': 'mm',
                          'soilmoist_p_tav': 'mm',                          
                          'tair_f_tavg': 'C',
                          'wind_f_tavg': 'm-s-1',
                          'tveg_tavg': 'mm'}
        elif step == 'daily':
            self.units = {'avgsurft_tavg': 'C',
                          'canopint_tavg': 'mm',
                          'evap_tavg': 'mm-day-1',
                          'lwdown_f_tavg': 'W-m-2',
                          'lwnet_tavg': 'W-m-2',
                          'psurf_f_tavg': 'kpa',
                          'qair_f_tavg': 'kg-kg',
                          'qg_tavg': 'W-m-2',
                          'qh_tavg': 'W-m-2',
                          'qle_tavg': 'W-m-2',
                          'qs_tavg': 'mm-day-1',
                          'qsb_tavg': 'mm-day-1',
                          'qsm_tavg': 'mm-day-1',
                          'rainf_f_tavg': 'mm-day-1',
                          'rainf_tavg': 'mm-day-1',
                          'swe_tavg': 'mm',
                          'swdown_f_tavg': 'W-m-2',
                          'swnet_tavg': 'W-m-2',
                          'snowf_tavg': 'mm',
                          'soilmoist_s_tav': 'mm',
                          'soilmoist_rz_ta': 'mm',
                          'soilmoist_p_tav': 'mm',                             
                          'tair_f_tavg': 'C',
                          'wind_f_tavg': 'm-s-1',
                          'tveg_tavg': 'mm'}
        elif step == 'monthly':
            self.units = {'avgsurft_tavg': 'C',
                          'canopint_tavg': 'mm',
                          'evap_tavg': 'mm-month-1',
                          'lwdown_f_tavg': 'W-m-2',
                          'lwnet_tavg': 'W-m-2',
                          'psurf_f_tavg': 'kpa',
                          'qair_f_tavg': 'kg-kg',
                          'qg_tavg': 'W-m-2',
                          'qh_tavg': 'W-m-2',
                          'qle_tavg': 'W-m-2',
                          'qs_tavg': 'mm-month-1',
                          'qsb_tavg': 'mm-month-1',
                          'qsm_tavg': 'mm-month-1',
                          'rainf_f_tavg': 'mm-month-1',
                          'rainf_tavg': 'mm-month-1',                          
                          'swe_tavg': 'mm',
                          'swdown_f_tavg': 'W-m-2',
                          'swnet_tavg': 'W-m-2',
                          'snowf_tavg': 'mm',
                          'soilmoist_s_tav': 'mm',
                          'soilmoist_rz_ta': 'mm',
                          'soilmoist_p_tav': 'mm',                             
                          'tair_f_tavg': 'C',
                          'wind_f_tavg': 'm-s-1',
                          'tveg_tavg': 'mm'}
        else:
            raise KeyError("The input time step is not supported")













https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_CLSM025_D.2.0%2F1948%2F01%2FGLDAS_CLSM025_D.A19480105.020.nc4&FORMAT=bmM0Lw&BBOX=37.969%2C-7.031%2C61.172%2C19.688&LABEL=GLDAS_CLSM025_D.A19480105.020.nc4.SUB.nc4&SHORTNAME=GLDAS_CLSM025_D&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.0
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_M.2.1%2F2001%2FGLDAS_NOAH025_M.A200101.021.nc4&FORMAT=Y29nLw&BBOX=25.313%2C-121.641%2C59.063%2C-82.969&LABEL=GLDAS_NOAH025_M.A200101.021.nc4.SUB.nc4&SHORTNAME=GLDAS_NOAH025_M&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Psurf_f_inst%2CWind_f_inst%2CTair_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_M.2.1%2F2000%2FGLDAS_NOAH025_M.A200001.021.nc4&FORMAT=bmM0Lw&BBOX=19.688%2C-119.531%2C50.625%2C-86.484&LABEL=GLDAS_NOAH025_M.A200001.021.nc4.SUB.nc4&SHORTNAME=GLDAS_NOAH025_M&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Psurf_f_inst%2CTair_f_inst%2CWind_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_M.2.1%2F2001%2FGLDAS_NOAH025_M.A200101.021.nc4&FORMAT=Y29nLw&BBOX=25.313%2C-121.641%2C59.063%2C-82.969&LABEL=GLDAS_NOAH025_M.A200101.021.nc4.SUB.nc4&SHORTNAME=GLDAS_NOAH025_M&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Psurf_f_inst%2CTair_f_inst%2CWind_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_3H.2.1%2F2000%2F001%2FGLDAS_NOAH025_3H.A20000101.0300.021.nc4&FORMAT=bmM0Lw&BBOX=-60%2C-180%2C90%2C180&LABEL=GLDAS_NOAH025_3H.A20000101.0300.021.nc4.SUB.nc4&SHORTNAME=GLDAS_NOAH025_3H&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Albedo_inst%2CLWdown_f_tavg%2CWind_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_M.2.1%2F2000%2FGLDAS_NOAH025_M.A200001.021.nc4&FORMAT=Y29nLw&BBOX=25.313%2C-121.641%2C59.063%2C-82.969&LABEL=GLDAS_NOAH025_M.A200001.021.nc4.SUB.tif&SHORTNAME=GLDAS_NOAH025_M&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Wind_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_M.2.1%2F2003%2FGLDAS_NOAH025_M.A200303.021.nc4&FORMAT=Y29nLw&BBOX=25.313%2C-121.641%2C59.063%2C-82.969&LABEL=GLDAS_NOAH025_M.A200303.021.nc4.SUB.tif&SHORTNAME=GLDAS_NOAH025_M&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Wind_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_M.2.0%2F2001%2FGLDAS_NOAH025_M.A200101.020.nc4&FORMAT=Y29nLw&BBOX=25.313%2C-121.641%2C59.063%2C-82.969&LABEL=GLDAS_NOAH025_M.A200101.020.nc4.SUB.tif&SHORTNAME=GLDAS_NOAH025_M&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.0&VARIABLES=Wind_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_M.2.1%2F2001%2FGLDAS_NOAH025_M.A200101.021.nc4&FORMAT=Y29nLw&BBOX=25.313%2C-121.641%2C59.063%2C-82.969&LABEL=GLDAS_NOAH025_M.A200101.021.nc4.SUB.tif&SHORTNAME=GLDAS_NOAH025_M&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Wind_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_3H.2.1%2F2000%2F001%2FGLDAS_NOAH025_3H.A20000101.1200.021.nc4&FORMAT=bmM0Lw&BBOX=-60%2C-180%2C90%2C180&LABEL=GLDAS_NOAH025_3H.A20000101.1200.021.nc4.SUB.nc4&SHORTNAME=GLDAS_NOAH025_3H&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Wind_f_inst

https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_3H.2.0%2F1948%2F001%2FGLDAS_NOAH025_3H.A19480101.2100.020.nc4&FORMAT=Y29nLw&BBOX=30.938%2C-104.062%2C51.328%2C-71.016&LABEL=GLDAS_NOAH025_3H.A19480101.2100.020.nc4.SUB.tif&SHORTNAME=GLDAS_NOAH025_3H&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.0&VARIABLES=Wind_f_inst
https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FGLDAS%2FGLDAS_NOAH025_3H.2.1%2F2000%2F001%2FGLDAS_NOAH025_3H.A20000101.0900.021.nc4&FORMAT=dGlmLw&BBOX=-60%2C-180%2C90%2C180&LABEL=GLDAS_NOAH025_3H.A20000101.0900.021.nc4.SUB.tif&SHORTNAME=GLDAS_NOAH025_3H&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=2.1&VARIABLES=Wind_f_inst