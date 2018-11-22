# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/FEWS

Description:
Use the FEWS.ETpot_daily function to
download and create monthly FEWS images in Gtiff format.
Data is available between 2003-01-01 till now. 

Examples:
from watertools.Collect import FEWS
FEWS.ETpot_daily(Dir='C:/Temp/', Startdate='2008-12-01', Enddate='2011-01-20',
           latlim=[-10, 30], lonlim=[-20, -10])
"""

from .ETpot_daily import main as ETpot_daily
__all__ = ['ETpot_daily']

__version__ = '0.1'
