# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MYD11

Description:
This module downloads MYD11 LST data from
http://e4ftl01.cr.usgs.gov/. Use the MYD11.LST_daily function to
download and create daily LST images in Gtiff format.
The data is available between 2000-02-18 till present.

Examples:
from watertools.Collect import MYD11
MYD11.LST_daily(Dir='C:/Temp3/', Startdate='2003-12-01', Enddate='2003-12-30',
           latlim=[41, 45], lonlim=[-8, -5])
"""

from .LST_daily import main as LST_daily
from .LST_8daily import main as LST_8daily

__all__ = ['LST_daily', 'LST_8daily']

__version__ = '0.1'
