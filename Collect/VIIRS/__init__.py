# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/VIIRS

Description:
This module downloads VIIRS LST data from

Examples:
from watertools.Collect import VIIRS
VIIRS.LST_daily(Dir='C:/Temp3/', Startdate='2003-12-01', Enddate='2003-12-30',
           latlim=[41, 45], lonlim=[-8, -5])
"""

from .LST_daily import main as LST_daily

__all__ = ['LST_daily']

__version__ = '0.1'
