# -*- coding: utf-8 -*-
'''
Authors: Tim Hessels
Module: Products/ETref
'''

# import watertools modules
from watertools.Collect import DEM
from watertools.Collect import CFSR
from watertools.Collect import GLDAS

def CollectData(Dir, Startdate, Enddate, latlim, lonlim, cores, CFSR_USE, GLDAS_USE, LANDSAF_USE):
    """
    This function Collect all the data needed for the ETref, by using wa.Collect functions.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax] (values must be between -60 and 60)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    cores -- amount of cores used for collecting data
    CFSR -- if 1 than CFSR data will be used
    """
    
    # Make latlim and lonlim a bit larger for GLDAS
    latlimGLDAS = [latlim[0] - 0.25, latlim[1] + 0.25]
    lonlimGLDAS = [lonlim[0] - 0.25, lonlim[1] + 0.25]
    
    # download DEM map (will only be done when it is not there yet)
    DEM.HydroSHED(Dir, latlimGLDAS, lonlimGLDAS, Waitbar = 0)
    
    if CFSR_USE == 1:
        # download CFSR data (will only be done when it is not there yet)
        CFSR.daily(Dir=Dir, Vars = ['dlwsfc','dswsfc','ulwsfc'], Startdate=Startdate, Enddate=Enddate, latlim=latlim, lonlim=lonlim, Waitbar = 0)

    if GLDAS_USE == 1:
        GLDAS.daily(Dir=Dir, Vars = ['swnet_tavg','lwnet_tavg'], Startdate=Startdate, Enddate=Enddate, latlim=latlimGLDAS, lonlim=lonlimGLDAS, SumMean=1, Min=0, Max=0, Waitbar = 0)

    # download GLDAS data (will only be done when it is not there yet)
    GLDAS.daily(Dir=Dir, Vars = ['tair_f_inst'], Startdate=Startdate, Enddate=Enddate, latlim=latlimGLDAS, lonlim=lonlimGLDAS, SumMean=0, Min=1, Max=1, Waitbar = 0)
    GLDAS.daily(Dir=Dir, Vars = ['psurf_f_inst','wind_f_inst','qair_f_inst'], Startdate=Startdate, Enddate=Enddate, latlim=latlimGLDAS, lonlim=lonlimGLDAS, SumMean=1, Min=0, Max=0, Waitbar = 0)

    return











