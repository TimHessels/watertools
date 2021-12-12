# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/MSGCPP
"""

import sys
from watertools.Collect.MSGCPP.DataAccess import DownloadData

def main(Dir, Startdate, Enddate, latlim, lonlim, Time = '', GMT_Offset = 0, Waitbar = 1):
    """
    This function downloads MOD11 daily data for the specified time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Time -- "12:45" Define the time, default is '' all 15min will be downloaded
    Waitbar -- 1 (Default) will print a waitbar
    """

    print('\nDownload 15min Clouds from MSGCPP for period %s till %s' %(Startdate, Enddate))    
        
    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Time, GMT_Offset, Waitbar, Type = "Cloud")

if __name__ == '__main__':
    main(sys.argv)