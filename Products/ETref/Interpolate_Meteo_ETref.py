# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Products/ETref
"""

# import general python modules
import os
import numpy as np


# import watertools modules
from watertools.General import data_conversions as DC
from watertools.General import raster_conversions as RC
from watertools.Products.ETref.SlopeInfluence_ETref import SlopeInfluence

def process_GLDAS(Tmax, Tmin, humidity, surface_pressure):
    """
    This function calculates the actual and saturated vapour pressure by using GLDAS data

    Keyword arguments:
    Tmax -- [] numpy array with the maximal temperature
    Tmin -- [] numpy array with the minimal temperature
    humidity -- [] numpy array with the humidity
    surface_pressure -- [] numpy array with the surface pressure
    """
	# calculate the average temparature based on FAO
    T_mean = (Tmax + Tmin) / 2

	# calculate the slope of saturation
    delta = 4098 * (0.6108  * np.exp(17.27 * T_mean / (T_mean + 237.3))) / (T_mean + 237.3)**2

    # calculate the saturation vapour pressure by using the max and min temparature
    e0_max = 0.6108 * np.exp(17.27 * Tmax / (Tmax + 237.3))
    e0_min = 0.6108 * np.exp(17.27 * Tmin / (Tmin + 237.3))

    # calculate the saturated vapour pressure
    es = (e0_max + e0_min) / 2

    # calculate the max and min relative humidity
    RH_max = np.minimum((1/0.622) * humidity * surface_pressure / e0_min, 1)*100
    RH_min = np.minimum((1/0.622) * humidity * surface_pressure / e0_max, 1)*100

    # calculate the actual vapour pressure
    ea = (e0_min * RH_max/100 + e0_max * RH_min/100) / 2

    return ea, es, delta

def lapse_rate(Dir,temperature_map, DEMmap):
    """
    This function downscales the GLDAS temperature map by using the DEM map

    Keyword arguments:
    temperature_map -- 'C:/' path to the temperature map
    DEMmap -- 'C:/' path to the DEM map
    """

    # calculate average altitudes corresponding to T resolution
    DEM_ave_out_name = os.path.join(Dir,'HydroSHED', 'DEM','DEM_ave.tif')

    if not os.path.exists(DEM_ave_out_name):
        dest = RC.reproject_dataset_example(DEMmap, temperature_map, method = 4)
        geo_out, proj, size_X, size_Y = RC.Open_array_info(temperature_map)
        DEM_ave_data = dest.GetRasterBand(1).ReadAsArray()
        DEM_ave_data[np.isnan(DEM_ave_data)] = 0
        DEM_ave_data = RC.gap_filling(DEM_ave_data, 0)
        DC.Save_as_tiff(DEM_ave_out_name, DEM_ave_data, geo_out, proj)

    else:
        DEM_ave_data = RC.Open_tiff_array(DEM_ave_out_name)
          
    dest = None

    # determine lapse-rate [degress Celcius per meter]
    lapse_rate_number = 0.0065

    # open maps as numpy arrays
    dest = RC.reproject_dataset_example(DEM_ave_out_name, DEMmap, method = 2)
    dem_avg=dest.GetRasterBand(1).ReadAsArray()
    dem_avg[dem_avg<0]=0
    dest = None

    # Open the temperature dataset
    dest = RC.reproject_dataset_example(temperature_map, DEMmap, method = 2)
    T=dest.GetRasterBand(1).ReadAsArray()
    T[np.isnan(T)] = -9999
    T = RC.gap_filling(T, -9999)

    
    dest = None

    # Open Demmap
    demmap = RC.Open_tiff_array(DEMmap)
    dem_avg[demmap<=0]=0
    demmap[demmap==-32768]=np.nan

    # calculate first part
    T = T + ((dem_avg-demmap) * lapse_rate_number)

    return T

def adjust_P(Dir, pressure_map, DEMmap):
    """
    This function downscales the GLDAS air pressure map by using the DEM map

    Keyword arguments:
    pressure_map -- 'C:/' path to the pressure map
    DEMmap -- 'C:/' path to the DEM map
    """

    # calculate average latitudes
    destDEMave = RC.reproject_dataset_example(DEMmap, pressure_map, method = 4)
    DEM_ave_out_name = os.path.join(Dir, 'HydroSHED', 'DEM','DEM_ave.tif')
    geo_out, proj, size_X, size_Y = RC.Open_array_info(pressure_map)
    DEM_ave_data = destDEMave.GetRasterBand(1).ReadAsArray()
    DC.Save_as_tiff(DEM_ave_out_name, DEM_ave_data, geo_out, proj)

    # open maps as numpy arrays
    dest = RC.reproject_dataset_example(DEM_ave_out_name, DEMmap, method = 2)
    dem_avg=dest.GetRasterBand(1).ReadAsArray()
    dest = None

    # open maps as numpy arrays
    dest = RC.reproject_dataset_example(pressure_map, DEMmap, method = 2)
    P=dest.GetRasterBand(1).ReadAsArray()
    dest = None

    demmap = RC.Open_tiff_array(DEMmap)
    dem_avg[demmap<=0]=0
    demmap[demmap==-32768]=np.nan

    # calculate second part
    P = P + (101.3*((293-0.0065*(demmap-dem_avg))/293)**5.26 - 101.3)

    os.remove(DEM_ave_out_name)

    return P

def slope_correct(down_short_hor, pressure, ea, DEMmap, DOY):
    """
    This function downscales the CFSR solar radiation by using the DEM map
    The Slope correction is based on Allen et al. (2006)
    'Analytical integrated functions for daily solar radiation on slope'

    Keyword arguments:
    down_short_hor -- numpy array with the horizontal downwards shortwave radiation
    pressure -- numpy array with the air pressure
    ea -- numpy array with the actual vapour pressure
    DEMmap -- 'C:/' path to the DEM map
    DOY -- day of the year
    """

    # Get Geo Info
    GeoT, Projection, xsize, ysize = RC.Open_array_info(DEMmap)

    minx = GeoT[0]
    miny = GeoT[3] + xsize*GeoT[4] + ysize*GeoT[5]

    x = np.flipud(np.arange(xsize)*GeoT[1] + minx + GeoT[1]/2)
    y = np.flipud(np.arange(ysize)*-GeoT[5] + miny + -GeoT[5]/2)

    # Calculate Extraterrestrial Solar Radiation [W m-2]
    demmap = RC.Open_tiff_array(DEMmap)
    demmap[demmap<0]=0

	# apply the slope correction
    Ra_hor, Ra_slp, sinb, sinb_hor, fi, slope, ID = SlopeInfluence(demmap,y,x,DOY)

    # Calculate atmospheric transmissivity
    Rs_hor = down_short_hor

    # EQ 39
    tau = Rs_hor/Ra_hor

    #EQ 41
    KB_hor = np.zeros(tau.shape) * np.nan

    indice = np.where(tau.flat >= 0.42)
    KB_hor.flat[indice] = 1.56*tau.flat[indice] -0.55

    indice = np.logical_and(tau.flat > 0.175, tau.flat < 0.42)
    KB_hor.flat[indice] = 0.022 - 0.280*tau.flat[indice] + 0.828*tau.flat[indice]**2 + 0.765*tau.flat[indice]**3

    indice = np.where(tau.flat <= 0.175)
    KB_hor.flat[indice] = 0.016*tau.flat[indice]

    # EQ 42
    KD_hor = tau - KB_hor

    Kt=0.7

    #EQ 18
    W = 0.14*ea*pressure + 2.1

    KB0 = 0.98*np.exp((-0.00146*pressure/Kt/sinb)-0.075*(W/sinb)**0.4)
    KB0_hor = 0.98*np.exp((-0.00146*pressure/Kt/sinb_hor)-0.075*(W/sinb_hor)**0.4)

    #EQ 34
    fB = KB0/KB0_hor * Ra_slp/Ra_hor
    fia = (1-KB_hor) * (1 + (KB_hor/(KB_hor+KD_hor))**0.5 * np.sin(slope/2)**3)*fi + fB*KB_hor

    Rs = Rs_hor*(fB*(KB_hor/tau) + fia*(KD_hor/tau) + 0.23*(1-fi))

    Rs[np.isnan(Rs)] = Rs_hor[np.isnan(Rs)]

    Rs_equiv = Rs / np.cos(slope)

    bias = np.nansum(Rs_hor)/np.nansum(Rs_equiv)

    return Rs_equiv, tau, bias
