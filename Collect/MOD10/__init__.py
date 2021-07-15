# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MOD10

Description:
This module downloads MOD10 SnowMask data from
https://n5eil01u.ecs.nsidc.org/MOST/MOD10A2.006. Use the MOD10.SnowMask_8daily function to
download and create 8 daily SnowMask images in Gtiff format.
The data is available between 2000-02-18 till present.

Examples:
from watertools.Collect import MOD10
MOD10.SnowMask_8daily(Dir='C:/Temp4/',Startdate='2003-12-01',Enddate='2003-12-20',
                                 latlim=[30,35],lonlim=[70,75])
"""

from .SnowMask_8daily import main as SnowMask_8daily
from .SnowMask_daily import main as SnowMask_daily

__all__ = ['SnowMask_8daily', 'SnowMask_daily']

__version__ = '0.1'
