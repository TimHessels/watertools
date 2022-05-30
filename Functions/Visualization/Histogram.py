# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:44:51 2020

@author: timhe
"""
import matplotlib.pyplot as plt
import numpy as np

def Show_array(Array, bins = 'auto'):
    
    rng = Array.flatten()
    rng[rng==-9999] = np.nan
    rng = rng[~np.isnan(rng)]
    plt.hist(rng, bins=bins)  # arguments are passed to np.histogram
    plt.title("Histogram with  %s bins" %bins)
    plt.show()
    return()