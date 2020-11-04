# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:09:38 2017

@author: tih
"""

import os
import watertools.General.raster_conversions as RC
import watertools.General.data_conversions as DC
import numpy as np
import pandas as pd

def Calc_Humidity(Temp_format, P_format, Hum_format, output_format, Startdate, Enddate, freq ="D"):

    folder_dir_out = os.path.dirname(output_format)
    
    if not os.path.exists(folder_dir_out):
        os.makedirs(folder_dir_out)

    Dates = pd.date_range(Startdate, Enddate, freq = freq)
    
    for Date in Dates:
        
        try:
            print(Date)
            
            Day = Date.day
            Month = Date.month
            Year = Date.year
            Hour = Date.hour
            Minute = Date.minute
            
            Tempfile_one = Temp_format.format(yyyy = Year, mm = Month, dd = Day, HH = Hour, MM = Minute)
            Presfile_one = P_format.format(yyyy = Year, mm = Month, dd = Day, HH = Hour, MM = Minute)
            Humfile_one = Hum_format.format(yyyy = Year, mm = Month, dd = Day, HH = Hour, MM = Minute)
            out_folder_one = output_format.format(yyyy = Year, mm = Month, dd = Day, HH = Hour, MM = Minute)
        
            geo_out, proj, size_X, size_Y = RC.Open_array_info(Tempfile_one)
            Tdata = RC.Open_tiff_array(Tempfile_one)
            Tdata[Tdata<-900]=-9999
            Pdata = RC.Open_tiff_array(Presfile_one)
            Hdata = RC.Open_tiff_array(Humfile_one)
            Pdata[Pdata<0]=-9999
            Hdata[Hdata<0]=-9999
            
            # gapfilling
            Tdata = RC.gap_filling(Tdata,-9999)
            Pdata = RC.gap_filling(Pdata,-9999)
            Hdata = RC.gap_filling(Hdata,-9999)
            
            if "_K_" in Temp_format:
                Tdata = Tdata - 273.15

            Esdata = 0.6108*np.exp((17.27*Tdata)/(Tdata+237.3))
            HumData = np.minimum((1.6077717*Hdata*Pdata/Esdata),1)*100
            HumData = HumData.clip(0,100)
                            
            DC.Save_as_tiff(out_folder_one,HumData,geo_out,"WGS84") 
            print("FINISHED relative humidity for %s" %Date)
        except:
            print("NO relative humidity for %s" %Date)

    return()   