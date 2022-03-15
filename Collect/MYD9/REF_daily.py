import sys
from watertools.Collect.MYD9.DataAccess import DownloadData


def main(Dir, Startdate, Enddate, latlim, lonlim, band = 1, resolution = "500m", cores=False, Waitbar = 1, hdf_library = None, remove_hdf = 1):
    """
    This function downloads MOD9 reflectance daily data for the specified time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
	 cores -- amount of cores used
    Waitbar -- 1 (Default) will print a waitbar
    """
    print('\nDownload daily MODIS Reflectance data for period %s till %s' %(Startdate, Enddate))
    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar, band, resolution, cores, hdf_library, remove_hdf)

if __name__ == '__main__':
    main(sys.argv)