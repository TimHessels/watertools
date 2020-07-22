# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 12:29:05 2020

@author: timhe
"""
import gdal
import numpy as np
import matplotlib.pyplot as plt

def Show_tif(image_file, Limits = None, Color = None):
    """
    This function plot a tiff array in the console

    Keyword arguments:
    image_file -- string
        Filename to the file that must be shown
    Limits -- [min, max] (Default = min and max image)
        User can define the limits of the colorbar
    Color -- string (Default = "viridis")
        User can define the wanted colormap, all options are listed here:
        https://matplotlib.org/examples/color/colormaps_reference.html
    """    
    dest = gdal.Open(image_file)
    Array = dest.GetRasterBand(1).ReadAsArray()
    Array[Array==-9999] = np.nan
    if Limits == None:
        Limits = [np.nanmin(Array), np.nanmax(Array)]
    
    if Color == None:
        Color = "viridis"
    
    plt.imshow(Array, cmap = Color, vmin=Limits[0], vmax=Limits[1])
    plt.colorbar()
    plt.show()
    
    return()
