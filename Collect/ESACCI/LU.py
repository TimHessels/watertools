# -*- coding: utf-8 -*-
import sys
import os
from watertools.Collect.ESACCI.DataAccess import DownloadData

def main(Dir, latlim, lonlim, Waitbar = 1):
    """
    This function downloads ESACCI daily data for a given variable, time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Waitbar -- 1 (Default) Will print a waitbar
    """
	
    # Create directory if not exists for the output
    output_folder = os.path.join(Dir, 'ESACCI', 'LU')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define the output map and create this if not exists
    nameEnd = os.path.join(output_folder, 'LU_ESACCI.tif')

    if not os.path.exists(nameEnd):

        # Create Waitbar
        if Waitbar == 1:
            print('\nDownload ESACCI landuse map')
            import watertools.Functions.Random.WaitbarConsole as WaitbarConsole
            total_amount = 1
            amount = 0
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

        # Download and process the data
        DownloadData(output_folder, latlim, lonlim, Waitbar)

        if Waitbar == 1:
            amount = 1
            WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    else:
        if Waitbar == 1:
            print("\nESACCI LU map already exists in output folder")

if __name__ == '__main__':
    main(sys.argv)

