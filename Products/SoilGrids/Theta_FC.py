# -*- coding: utf-8 -*-
'''
Authors: Tim Hessels
Module: Products/SoilGrids
'''

# import general python modules
import os
import numpy as np

# import WA+ modules
from watertools.General import data_conversions as DC
from watertools.General import raster_conversions as RC

def Topsoil(Dir, latlim, lonlim):
    """
    This function calculates the Theta Field Capacity soil characteristic (15cm)

    Keyword arguments:
    Dir -- 'C:/' path to the WA map
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    """
	
    print('/nCreate Theta Field Capacity map of the topsoil from SoilGrids')
    
    # Define parameters to define the topsoil
    SL = "sl3"
	
    Calc_Property(Dir, latlim, lonlim, SL)
	
    return

def Subsoil(Dir, latlim, lonlim):
    """
    This function calculates the Theta Field Capacity characteristic (100cm)

    Keyword arguments:
    Dir -- 'C:/' path to the WA map
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    """
	
    print('/nCreate Theta Field Capacity map of the subsoil from SoilGrids')
    
	 # Define parameters to define the subsoil	
    SL = "sl6"
	
    Calc_Property(Dir, latlim, lonlim, SL)
	
    return
	
def Calc_Property(Dir, latlim, lonlim, SL):	

    import watertools
    
    # Define level
    if SL == "sl3":
        level = "Topsoil"
    elif SL == "sl6":
        level = "Subsoil" 
    
    # check if you need to download
    filename_out_thetasat = os.path.join(Dir, 'SoilGrids', 'Theta_Sat' ,'Theta_Sat2_%s_SoilGrids_kg-kg.tif' %level)
    if not os.path.exists(filename_out_thetasat):
       if SL == "sl3":
           watertools.Products.SoilGrids.Theta_Sat2.Topsoil(Dir, latlim, lonlim)
       elif SL == "sl6":
           watertools.Products.SoilGrids.Theta_Sat2.Subsoil(Dir, latlim, lonlim)
           
    filename_out_thetares = os.path.join(Dir, 'SoilGrids', 'Theta_Res' ,'Theta_Res_%s_SoilGrids_kg-kg.tif' %level)
    if not os.path.exists(filename_out_thetares):
       if SL == "sl3":
           watertools.Products.SoilGrids.Theta_Res.Topsoil(Dir, latlim, lonlim)
       elif SL == "sl6":
           watertools.Products.SoilGrids.Theta_Res.Subsoil(Dir, latlim, lonlim)           
           
    filename_out_n_genuchten = os.path.join(Dir, 'SoilGrids', 'N_van_genuchten' ,'N_genuchten_%s_SoilGrids_-.tif' %level)
    if not os.path.exists(filename_out_n_genuchten):
       if SL == "sl3":
           watertools.Products.SoilGrids.n_van_genuchten.Topsoil(Dir, latlim, lonlim)
       elif SL == "sl6":
           watertools.Products.SoilGrids.n_van_genuchten.Subsoil(Dir, latlim, lonlim)              

    filedir_out_thetafc = os.path.join(Dir, 'SoilGrids', 'Theta_FC')
    if not os.path.exists(filedir_out_thetafc):
        os.makedirs(filedir_out_thetafc)   
              
    # Define theta field capacity output
    filename_out_thetafc = os.path.join(filedir_out_thetafc, 'Theta_FC2_%s_SoilGrids_cm3-cm3.tif' %level)

    if not os.path.exists(filename_out_thetafc):
            
        # Get info layer
        geo_out, proj, size_X, size_Y = RC.Open_array_info(filename_out_thetasat)
        
        # Open dataset
        theta_sat = RC.Open_tiff_array(filename_out_thetasat)
        theta_res = RC.Open_tiff_array(filename_out_thetares)
        n_genuchten = RC.Open_tiff_array(filename_out_n_genuchten)
        
        # Calculate theta field capacity
        theta_FC = np.ones(theta_sat.shape) * -9999   
        #theta_FC = np.where(theta_sat < 0.301, 0.042, np.arccosh(theta_sat + 0.7) - 0.32 * (theta_sat + 0.7) + 0.2)        
        #theta_FC = np.where(theta_sat < 0.301, 0.042, -2.95*theta_sat**2+3.96*theta_sat-0.871)   
        
        theta_FC = theta_res +  (theta_sat - theta_res)/(1+(0.02 * 200)**n_genuchten)**(1-1/n_genuchten)                
        # Save as tiff
        DC.Save_as_tiff(filename_out_thetafc, theta_FC, geo_out, proj)

    return           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           