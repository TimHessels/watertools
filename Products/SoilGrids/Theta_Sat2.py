# -*- coding: utf-8 -*-
'''
Authors: Tim Hessels
Module: Products/SoilGrids
'''

# import general python modules
import os
import gdal
import numpy as np

# import WA+ modules
from watertools.General import data_conversions as DC
from watertools.General import raster_conversions as RC

def Topsoil(Dir, latlim, lonlim, GF = True):
    """
    This function calculates the topsoil saturated soil characteristic (15cm)

    Keyword arguments:
    Dir -- 'C:/' path to the WA map
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    """
	
    print('/nCreate Theta Saturated map of the topsoil from SoilGrids')
    
    # Define parameters to define the topsoil
    SL = "sl3"
	
    Calc_Property(Dir, latlim, lonlim, SL, GF)
	
    return

def Subsoil(Dir, latlim, lonlim, GF = True):
    """
    This function calculates the subsoil saturated soil characteristic (100cm)

    Keyword arguments:
    Dir -- 'C:/' path to the WA map
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    """
	
    print('/nCreate Theta Saturated map of the subsoil from SoilGrids')
    
	# Define parameters to define the subsoil	
    SL = "sl6"
	
    Calc_Property(Dir, latlim, lonlim, SL, GF)
	
    return
	
def Calc_Property(Dir, latlim, lonlim, SL, GF = True):	
	
    import watertools.Collect.SoilGrids as SG

    # Download needed layers
    SG.Clay_Content(Dir, latlim, lonlim, level=SL)
    SG.Silt_Content(Dir, latlim, lonlim, level=SL)    
    SG.Organic_Carbon_Content(Dir, latlim, lonlim, level=SL)
    SG.Bulk_Density(Dir, latlim, lonlim, level=SL)
    
    # Define path to layers
    filename_clay = os.path.join(Dir, 'SoilGrids', 'Clay_Content' ,'ClayContentMassFraction_%s_SoilGrids_percentage.tif' %SL)
    filename_om = os.path.join(Dir, 'SoilGrids', 'Soil_Organic_Carbon_Content' ,'SoilOrganicCarbonContent_%s_SoilGrids_g_kg.tif' %SL)
    filename_bulkdensity = os.path.join(Dir, 'SoilGrids', 'Bulk_Density' ,'BulkDensity_%s_SoilGrids_kg-m-3.tif' %SL)
    filename_silt = os.path.join(Dir, 'SoilGrids', 'Silt_Content' ,'SiltContentMassFraction_%s_SoilGrids_percentage.tif' %SL)
    
    # Define path for output
    if SL == "sl3":
        level = "Topsoil"
    elif SL == "sl6":
        level = "Subsoil"        

    filedir_out_densbulk = os.path.join( Dir, 'SoilGrids', 'Bulk_Density')
    if not os.path.exists(filedir_out_densbulk):
        os.makedirs(filedir_out_densbulk)
    filedir_out_thetasat = os.path.join(Dir, 'SoilGrids', 'Theta_Sat')
    if not os.path.exists(filedir_out_thetasat):
        os.makedirs(filedir_out_thetasat)    
    
    #filename_out_densbulk = os.path.join(filedir_out_densbulk ,'Bulk_Density_%s_SoilGrids_g-cm-3.tif' %level)
    filename_out_thetasat = os.path.join(filedir_out_thetasat,'Theta_Sat2_%s_SoilGrids_kg-kg.tif' %level)

    #if not (os.path.exists(filename_out_densbulk) and os.path.exists(filename_out_thetasat)):
    if not os.path.exists(filename_out_thetasat):
        
        # Open datasets
        dest_clay = gdal.Open(filename_clay)
        dest_om = gdal.Open(filename_om)
        dest_bulk = gdal.Open(filename_bulkdensity)
        dest_silt = gdal.Open(filename_silt)

        # Open Array info
        geo_out, proj, size_X, size_Y = RC.Open_array_info(filename_clay)    
        
        # Open Arrays
        Clay = dest_clay.GetRasterBand(1).ReadAsArray()
        OM = dest_om.GetRasterBand(1).ReadAsArray()
        Silt = dest_silt.GetRasterBand(1).ReadAsArray()        

        Clay = np.float_(Clay)
        '''
        try:
           Clay = RC.gap_filling(Clay, 0, method = 1)
        except:
            pass
        '''        
        OM = np.float_(OM)
        '''
        try:
           OM = RC.gap_filling(OM, 0, method = 1)
        except:
            pass
        '''
        Silt = np.float_(Silt)
        '''        
        try:
           Silt = RC.gap_filling(Silt, 0, method = 1)
        except:
            pass
        '''        
        Clay[Clay>100]=np.nan                
        Silt[Silt>100]=np.nan
        OM = OM/1000
        
        Clay[Clay==0]=np.nan  
        Silt[Silt==0]=np.nan  
        OM[OM==0]=np.nan  
        
        # Calculate bulk density
        bulk_dens1 = dest_bulk.GetRasterBand(1).ReadAsArray()
        bulk_dens1 = bulk_dens1/1000 
        bulk_dens1 = np.float_(bulk_dens1)
        '''           
        try:
           bulk_dens1 = RC.gap_filling(bulk_dens1, 0, method = 1)
        except:
            pass
        '''                   
        bulk_dens2 = 1/(0.6117 + 0.3601 * Clay/100 + 0.002172 * np.power(OM * 100, 2)+ 0.01715 * np.log(OM * 100))
    
        '''           
        try:
           bulk_dens2 = RC.gap_filling(bulk_dens2, 0, method = 1)
        except:
            pass       
        '''           
        bulk_dens = np.where(bulk_dens1>bulk_dens2, bulk_dens1, bulk_dens2)
        
        '''
        # Oude methode gebaseerd op Schenost, Sinowski & Priesack (1996) 
        
        # Calculate theta saturated
        theta_sat = 0.85 * (1- (bulk_dens/2.65)) + 0.13 * Clay/100
        '''
        # Nieuwe methode gebaseerd op Toth et al (2014)
        
        # Calculate silt fraction based on clay fraction
        #Silt_fraction = 0.7 * (Clay/100) ** 2 + 0.308 * Clay/100
        Silt_fraction = Silt/100
        
        # Calculate theta sat
        theta_sat = 0.8308 - 0.28217 * bulk_dens + 0.02728 * Clay/100 + 0.0187 * Silt_fraction
        theta_sat[Clay==0] = -9999
        if GF == True:
            try:
               theta_sat = RC.gap_filling(theta_sat, -9999, method = 1)
            except:
                pass   

        # Save data
        #DC.Save_as_tiff(filename_out_densbulk, bulk_dens, geo_out, "WGS84")
        DC.Save_as_tiff(filename_out_thetasat, theta_sat, geo_out, "WGS84")
        
    return()    

