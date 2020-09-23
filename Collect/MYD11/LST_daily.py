# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MYD11
"""

import sys
from watertools.Collect.MYD11.DataAccess import DownloadData

def main(Dir, Startdate, Enddate, latlim, lonlim, Waitbar = 1, cores=False, hdf_library = None, remove_hdf = 1, angle_info = 0, time_info = 0):
    """
    This function downloads MYD11 daily data for the specified time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Waitbar -- 1 (Default) will print a waitbar
    hdf_library -- string, if all the hdf files are already stored on computer
                    define directory to the data here
    remove_hdf -- 1 (Default), if 1 remove all the downloaded hdf files in the end
    """

    print('\nDownload daily MODIS land surface temperature data for period %s till %s' %(Startdate, Enddate))
    TimeStep = 1
    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, TimeStep, Waitbar, cores, hdf_library, remove_hdf, angle_info, time_info)

if __name__ == '__main__':
    main(sys.argv)