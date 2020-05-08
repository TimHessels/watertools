# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Contact: timhessels@hotmail.com
Repository: https://github.com/TimHessels/watertools
Module: Collect/SoilGrids


Description:


Examples:



"""
from .Bulk_Density import main as Bulk_Density
from .Clay_Content import main as Clay_Content
from .Organic_Carbon_Content import main as Organic_Carbon_Content
from .Organic_Carbon_Density import main as Organic_Carbon_Density
from .Sand_Content import main as Sand_Content
from .Silt_Content import main as Silt_Content
from .Nitrogen import main as Nitrogen
from .Soil_pH import main as Soil_pH

__all__ = ['Bulk_Density','Clay_Content','Organic_Carbon_Content','Organic_Carbon_Density','Sand_Content','Silt_Content','Nitrogen','Soil_pH']

__version__ = '0.1'
