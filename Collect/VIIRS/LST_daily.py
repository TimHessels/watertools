# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 13:16:00 2022

@author: timhe
"""

import sys
from watertools.Collect.VIIRS.DataAccess import DownloadData

def main(Dir, Startdate, Enddate, latlim, lonlim, Waitbar = 1):
    """
    This function downloads daily LST data from VIIRS
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Waitbar -- 1 (Default) will print a waitbar
    """

    print('\nDownload daily VIIRS LST data for period %s till %s' %(Startdate, Enddate))
    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar)

if __name__ == '__main__':
    main(sys.argv)