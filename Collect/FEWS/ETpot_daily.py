# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 07:54:17 2017

@author: tih
"""
import sys
from watertools.Collect.SSEBop.DataAccess import DownloadData

def main(Dir, Startdate='', Enddate='', latlim=[-59.17, 80], lonlim=[-180, 180], version = "V6", Waitbar = 1):
    """
    This function downloads daily FEWS Potential Evapotranspiration data

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax] (values must be between -59.17 and 80)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    """
    print('\nDownload daily FEWS potential evapotranspiration data for the period %s till %s' %(Startdate, Enddate))

    # Download data
    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar, version, "daily", Product = "ETpot")
                 
if __name__ == '__main__':
    main(sys.argv)
