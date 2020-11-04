# -*- coding: utf-8 -*-
import sys
from watertools.Collect.GEOS.DataAccess import DownloadData

def main(Dir, Vars, Startdate, Enddate, latlim, lonlim, Periods = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24], Waitbar = 1):
    """
    This function downloads GEOS inst data for a given variable, time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Vars -- [''] 
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
	Periods -- [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24] Period that needs to be downloaded. 1 period is 1 hour starting from noon
    Waitbar -- 1 (Default) Will print a waitbar
    """
   
    for Var in Vars:

        for Period in Periods:
		
            if Waitbar == 1:
                print('\nDownloading hourly GEOS %s data for the period %s till %s, Period = %s' %(Var, Startdate, Enddate, Period))

            # Download data
            DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, "hourly", Period, Waitbar, data_type = ["mean"])

if __name__ == '__main__':
    main(sys.argv)
