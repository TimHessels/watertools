# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/NLDAS


Description:
This script automatically downloads NLDAS data with a 0.125 degree
resolution for different extends based on their opendap server. The time
interval available are: three-hourly ('3hour'), daily ('day'), and monthly
('month'). A list of the variables can be printed with the command:
    NLDAS.VarInfo('daily').descriptions.keys()
Futher information of the variable can be printed with the following commands:
    NLDAS.VarInfo('daily').descriptions['evap']
    NLDAS.VarInfo('daily').units['evap']
    NLDAS.VarInfo('daily').names['evap']
    NLDAS.VarInfo('daily').factors['evap']

Examples:
from watertools.Collect import NLDAS
NLDAS.FORA_hourly(Dir='C:/Temp/', Vars=['tmp2m','pressfc'], Startdate='2004-12-20', Enddate='2005-01-10',
                   latlim=[38, 41], lonlim=[-76, -73], Periods=[4, 5])
NLDAS.FORA_daily(Dir='C:/Temp/', Vars=['tmp2m'], Startdate='2004-12-20', Enddate='2005-01-01',
            latlim=[38, 41], lonlim=[-76, -73],
            SumMean=1, Min=1, Max=1)
NLDAS.FORA_monthly(Dir='C:/TempNLDAS', Vars=['tmp2m'], Startdate='2004-12-20', Enddate='2005-03-10',latlim=[38, 41], lonlim=[-76, -73])
"""

from .FORA_DataAccess import VariablesInfo as FORA_VarInfo
from .FORA_daily import main as FORA_daily
from .FORA_monthly import main as FORA_monthly
from .FORA_hourly import main as FORA_hourly
from .DataAccess import VariablesInfo as VarInfo
from .daily import main as daily
from .monthly import main as monthly
from .hourly import main as hourly

__all__ = ['FORA_VarInfo', 'FORA_daily', 'FORA_monthly', 'FORA_hourly', 'VarInfo', 'daily', 'monthly', 'hourly']

__version__ = '0.1'
