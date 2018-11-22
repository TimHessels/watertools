# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/JRC

Description:
This module downloads JRC water occurrence data from http://storage.googleapis.com/global-surface-water/downloads/.
Use the JRC.Occurrence function to
download and create a water occurrence image in Gtiff format.
The data represents the period 1984-2015.

Examples:
from watertools.Collect import JRC
JRC.Occurrence(Dir='C:/Temp3/', latlim=[41, 45], lonlim=[-8, -5])
"""

from .Occurrence import main as Occurrence

__all__ = ['Occurrence']

__version__ = '0.1'
