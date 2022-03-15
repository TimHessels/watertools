# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MYD9

Description:
This module downloads MYD9 reflectance data from
https://e4ftl01.cr.usgs.gov/. Use the MOD9.REF_daily function to
download and create daily Reflectance images in Gtiff format.
The data is available between 2000-02-18 till present.

Examples:
from watertools.Collect import MYD9
MYD9.REF_daily(Dir='C:/Temp3/', Startdate='2003-12-01', Enddate='2003-12-20',
           latlim=[41, 45], lonlim=[-8, -5])
"""

from .REF_daily import main as REF_daily

__all__ = ['REF_daily']

__version__ = '0.1'
