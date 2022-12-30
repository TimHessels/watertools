# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 12:29:05 2020

@author: timhe
"""
import os
from osgeo import gdal
import glob
import numpy as np
import matplotlib.pyplot as plt

def Show_tif(image_file, Limits = None, Color = None, Invert = False, output_format = None, title = None):
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

    directory = os.path.dirname(image_file)
    os.chdir(directory)
    file = glob.glob(image_file)[0]
    
    dest = gdal.Open(os.path.join(directory, file))
    Array = dest.GetRasterBand(1).ReadAsArray()
    
    Array = np.float_(Array)
    Array[Array==-9999] = np.nan
    if Limits == None:
        Limits = [np.nanmin(Array), np.nanmax(Array)]
    
    if Color == None:
        Color = "viridis"
    if Invert == True:
        Color = plt.cm.get_cmap(Color).reversed()
        
    plt.imshow(Array, cmap = Color, vmin=Limits[0], vmax=Limits[1])
    if title != None:
        plt.title(title)
        
    plt.colorbar()
    
    if output_format != None:
        directory_out = os.path.dirname(output_format)
        if not os.path.exists(directory_out):
            os.makedirs(directory_out)
        
        plt.savefig(output_format, bbox_inches='tight')
    else:
        plt.show()
        
    plt.close()
    
    return()


def Show_array(image_file, Limits = None, Color = None):
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

    Array = image_file
    Array = np.float_(Array)
    Array[Array==-9999] = np.nan
    if Limits == None:
        Limits = [np.nanmin(Array), np.nanmax(Array)]
    
    if Color == None:
        Color = "viridis"
    
    plt.imshow(Array, cmap = Color, vmin=Limits[0], vmax=Limits[1])
    plt.colorbar()
    plt.show()
    
    return()