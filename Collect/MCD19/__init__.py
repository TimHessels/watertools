# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MCD19

Description:
This module downloads MCD19 reflectance data from
https://e4ftl01.cr.usgs.gov/. Use the MCD19.Albedo_8daily function to
download and create 8-daily Reflectance images in Gtiff format.
The data is available between 2000-02-18 till present.

Examples:
from watertools.Collect import MCD19
MCD19.Albedo_8daily(Dir='C:/Temp3/', Startdate='2003-12-01', Enddate='2003-12-20',
           latlim=[41, 45], lonlim=[-8, -5])
"""

from .Albedo_8daily import main as Albedo_8daily
from .Albedo_daily import main as Albedo_daily

__all__ = ['Albedo_8daily', 'Albedo_daily']

__version__ = '0.1'
