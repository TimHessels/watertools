# -*- coding: utf-8 -*-
"""
    Authors: Tim Hessels
    Module: Products/GridSoils

    Description:
    This product will calculate soil properties by using the GridSoils as basis. 
	Different soil characteristics can be estimated.
	The formulas are taken from the SoMoi model and the SoilGrids are taken from:
	ftp.soilgrids.org
"""

from watertools.Products.SoilGrids import K_Sat
from watertools.Products.SoilGrids import Theta_FC
from watertools.Products.SoilGrids import Theta_Sat
from watertools.Products.SoilGrids import Theta_Sat2
from watertools.Products.SoilGrids import Theta_Res
from watertools.Products.SoilGrids import Water_Holding_Capacity
from watertools.Products.SoilGrids import n_van_genuchten

__all__ = ['K_Sat', 'Theta_FC', 'Theta_Sat', 'Theta_Sat2', 'Theta_Res', 'Water_Holding_Capacity', 'n_van_genuchten']

__version__ = '0.1'
