# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/DEM


Description:
This module downloads DEM data from http://earlywarning.usgs.gov/hydrodata/.
Use the DEM functions to download and create DEM images in Gtiff format.

Examples:
from watertools.Collect import DEM
DEM.HydroSHED(Dir='C:/TempDEM4/', latlim=[29, 32], lonlim=[-113, -109])
"""
from .HydroSHED import main as HydroSHED
from .HydroSHED_Dir import main as HydroSHED_Dir
from .HydroSHED_Acc import main as HydroSHED_Acc
from .SRTM import main as SRTM

__all__ = ['HydroSHED','HydroSHED_Dir', 'HydroSHED_Acc', 'SRTM']

__version__ = '0.1'
