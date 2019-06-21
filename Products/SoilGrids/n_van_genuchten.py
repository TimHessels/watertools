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
    This function calculates the n van genuchten soil characteristic (15cm)

    Keyword arguments:
    Dir -- 'C:/' path to the WA map
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    """
	
    print('/nCreate n van genuchten map of the topsoil from SoilGrids')
    
    # Define parameters to define the topsoil
    SL = "sl3"
	
    Calc_Property(Dir, latlim, lonlim, SL)
	
    return

def Subsoil(Dir, latlim, lonlim):
    """
    This function calculates the subsoil n van genuchten soil characteristic (100cm)

    Keyword arguments:
    Dir -- 'C:/' path to the WA map
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    """
	
    print('/nCreate n van genuchten map of the subsoil from SoilGrids')
    
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

    filedir_out_n_genuchten = os.path.join(Dir, 'SoilGrids', 'N_van_genuchten')
    if not os.path.exists(filedir_out_n_genuchten):
        os.makedirs(filedir_out_n_genuchten)   
             
    # Define n van genuchten output
    filename_out_ngenuchten = os.path.join(filedir_out_n_genuchten ,'N_genuchten_%s_SoilGrids_-.tif' %level)

    if not os.path.exists(filename_out_ngenuchten):
            
        # Get info layer
        geo_out, proj, size_X, size_Y = RC.Open_array_info(filename_out_thetasat)
        
        # Open dataset
        theta_sat = RC.Open_tiff_array(filename_out_thetasat)
        
        # Calculate n van genuchten
        n_van_genuchten = np.ones(theta_sat.shape) * -9999   
        n_van_genuchten = 166.63*theta_sat**4-387.72*theta_sat**3+340.55*theta_sat**2-133.07*theta_sat+20.739
		
        # Save as tiff
        DC.Save_as_tiff(filename_out_ngenuchten, n_van_genuchten, geo_out, proj)
    return           
           
           
	
	


