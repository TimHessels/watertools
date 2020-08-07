# -*- coding: utf-8 -*-
import sys
from watertools.Collect.NLDAS.FORA_DataAccess import DownloadData


def main(Dir, Vars, Startdate, Enddate, latlim, lonlim, cores=False,
         Periods=[1, 2, 3, 4, 5, 6, 7, 8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24], Waitbar = 1):
    """
    This function downloads NLDAS three-hourly data for a given variable, time
    interval, spatial extent, and day period.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Var --  ['tmp2m','spfh2m'] Variable code. Run: VariablesInfo('day').descriptions.keys()
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Periods -- List of numbers from 1 to 24 (e.g. [1,4,5,8]). Stands for the
               period of hour of a day as follows:
                    Period       Hours
                      1      00:00 - 01:00
                      2      01:00 - 02:00
                      3      02:00 - 03:00
                      4      03:00 - 04:00
                      etc
    Waitbar -- 1 (Default) Will print a waitbar

    """
    for Var in Vars:

        if Waitbar == 1:
            print('\nDownloading hourly NLDAS %s data for the period %s till %s' %(Var, Startdate, Enddate))

        # Download data
        DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, Waitbar, cores,
					 TimeCase='hourly', CaseParameters=Periods)    # Download data

if __name__ == '__main__':
    main(sys.argv)
