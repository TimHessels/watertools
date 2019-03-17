# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Mon Feb 11 08:57:38 2019
"""

import sys
from watertools.Collect.L30.DataAccess import DownloadData


def main(Dir, Startdate, Enddate, S2_tiles):
                """
                Keyword arguments:
                Dir -- 'C:/file/to/path/'
                Startdate -- 'yyyy-mm-dd'
                Enddate -- 'yyyy-mm-dd'
                S2_tiles -- Tiles that you want to download
                """
                print('\nDownload L30 Harmonized data for the period %s till %s' %(Startdate, Enddate))
                
                for S2_tile in S2_tiles:
                    # Download data
                    DownloadData(Dir, Startdate, Enddate, S2_tile)

if __name__ == '__main__':
    main(sys.argv)