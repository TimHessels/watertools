# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/TRMM


Description:
This module downloads TRMM3B42 V7 (daily) and TRMM3B43 V7 (monthly) data from
ftp://disc2.nascom.nasa.gov. Use the TRMM.daily or TRMM.monthly functions to
download and create daily or monthly TRMM images in Gtiff format.
The TRMM data is available since 1998-01-01 till 2015-04-30

Examples:
from watertools.Collect import TRMM
TRMM.daily(Dir='C:/Temp/', Startdate='1999-02-01', Enddate='1999-02-28',
           latlim=[-10, 30], lonlim=[-20, 120])
TRMM.monthly(Dir='C:/Temp/', Startdate='1999-02-01', Enddate='1999-02-28',
             latlim=[-10, 30], lonlim=[-20, 120])
"""

from watertools.Collect.RFE.daily import main as daily
from watertools.Collect.RFE.monthly import main as monthly

__all__ = ['daily', 'monthly']

__version__ = '0.1'
