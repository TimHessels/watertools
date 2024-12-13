# -*- coding: utf-8 -*-
'''
Authors: Tim Hessels
Module: Products/ETref
'''
# import general python modules
import sys
import numpy as np

# import watertools modules
from watertools.Products.ETref.CollectDataETref import CollectData
from watertools.Products.ETref.CollectLANDSAFETref import CollectLANDSAF
from watertools.Products.ETref.SetVarETref import SetVariables

def main(Dir, Startdate = '', Enddate = '',
         latlim = [-60, 60], lonlim = [-180, 180], pixel_size = False, cores = False, LANDSAF_USE =  0, CFSR_USE =  0, GLDAS_USE =  1, SourceLANDSAF=  '', Waitbar = 1):
    """
    This function creates ETref (daily) data based on Hydroshed, GLDAS, and (CFSR/LANDSAF)

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax] (values must be between -60 and 60)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    pixel_size -- The output pixel size
    cores -- The number of cores used to run the routine.
             It can be 'False' to avoid using parallel computing
             routines.
    LANDSAF -- if LANDSAF data must be used it is 1
    SourceLANDSAF -- the path to the LANDSAF files
    Waitbar -- 1 (Default) will print the waitbar
    """
    if Waitbar == 1:
        print('Create daily Reference ET data for the period %s till %s' %(Startdate, Enddate))

    # Correct latitude and longitude if needed
    if latlim[0] < -60 or latlim[1] > 60:
        print('Latitude above 60N or below 60S is not possible.'
                        ' Value set to maximum')
        latlim[0] = np.max(latlim[0], -60)
        latlim[1] = np.min(lonlim[1], 60)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W.'
                       ' Now value is set to maximum')
        lonlim[0] = np.max(latlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)

    # Download data (using the wa.Collect scripts)
    CollectData(Dir, Startdate, Enddate, latlim, lonlim, cores, CFSR_USE, GLDAS_USE, LANDSAF_USE)

    # Process LANDSAF data if needed
    if LANDSAF_USE == 1:
        CollectLANDSAF(SourceLANDSAF, Dir, Startdate, Enddate, latlim, lonlim)

    # Set up the variables and calculates ETref hereafter
    SetVariables(Dir, Startdate, Enddate, latlim, lonlim, pixel_size, cores, LANDSAF_USE, CFSR_USE, GLDAS_USE, Waitbar)


# if __name__ == '__main__':
#     main(sys.argv)
