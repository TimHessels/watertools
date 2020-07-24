# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Products


Description:
This module contains scripts used to create WA+ products (data directly from web).

Products                      Dates                             Password                        
ETens (monthly)               2003/01/01-2014/12/31             WA+ FTP
ETref (daily)                 2000/01/01-now                    NASA
ETref (monthly)               2000/01/01-now                    NASA

Examples:
from watertools import Products
help(Products)
dir(Products)
"""

from watertools.Products import ETref
from watertools.Products import ETens
from watertools.Products import SoilGrids
from watertools.Products import RH

__all__ = ['ETref', 'ETens', 'SoilGrids', 'RH']

__version__ = '0.1'