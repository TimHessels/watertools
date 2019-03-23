# -*- coding: utf-8 -*-

# General modules
import numpy as np
import calendar
import os
import pandas as pd
import requests
from joblib import Parallel, delayed

# Water Accounting modules
from watertools import WebAccounts
import watertools.General.data_conversions as DC

def DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, Waitbar, cores,
                 TimeCase, CaseParameters):
    """
    This function downloads NLDAS Forcing data hourly, daily or monthly data

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Var -- 'wind_f_inst' : (string) For all variable codes: VariablesInfo('day').descriptions.keys()
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    cores -- 1....8
    CaseParameters -- See files: three_hourly.py, daily.py, and monthly.py
    """

    # Load factors / unit / type of variables / accounts
    VarInfo = VariablesInfo(TimeCase)
    username, password = WebAccounts.Accounts(Type = 'NASA')

    # Set required data for the three hourly option
    if TimeCase == 'hourly':

        # Define output folder and create this one if not exists
        path = os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS',
                            TimeCase, Var)

        if not os.path.exists(path):
            os.makedirs(path)

        # Startdate if not defined
        sd_date = '1979-01-02'

        # Define Time frequency
        TimeFreq = 'D'

        # Define URL by using personal account
        #url = 'http://%s:%s@hydro1.gesdisc.eosdis.nasa.gov:80/dods/GLDAS_NOAH025SUBP_3H' %(username,password)
        url = 'https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_NOAH0125_H.002' #%(username,password)

        # Name the definition that will be used to obtain the data
        RetrieveData_fcn = RetrieveData_hourly

    # Set required data for the daily option
    elif TimeCase == 'daily':

        # seperate the daily case parameters
        SumMean, Min, Max = CaseParameters

        # Define output folder and create this one if not exists
        path = {'mean': os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS',
                                     TimeCase, Var, 'mean'),
                'min': os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS',
                                    TimeCase, Var, 'min'),
                'max': os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS',
                                    TimeCase, Var, 'max')}
        selected = np.array([SumMean, Min, Max])
        types = np.array(('mean', 'min', 'max'))[selected == 1]
        CaseParameters = [selected, types]
        for i in range(len(types)):
            if not os.path.exists(path[types[i]]):
                os.makedirs(path[types[i]])

        # Startdate if not defined
        sd_date = '1979-01-02'

        # Define Time frequency
        TimeFreq = 'D'

        # Define URL by using personal account
        #url = 'http://%s:%s@hydro1.gesdisc.eosdis.nasa.gov:80/dods/GLDAS_NOAH025SUBP_3H' %(username,password)
        url = 'https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_NOAH0125_H.002' #%(username,password)

        # Name the definition that will be used to obtain the data
        RetrieveData_fcn = RetrieveData_daily

    # Set required data for the monthly option
    elif TimeCase == 'monthly':

        # Define output folder and create this one if not exists
        path = os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS',
                            TimeCase, Var)
        if not os.path.exists(path):
            os.makedirs(path)
        CaseParameters = []

        # Startdate if not defined
        sd_date = '1979-02-01'

        # Define Time frequency
        TimeFreq = 'MS'

        # Define URL by using personal account
        #url = 'http://%s:%s@hydro1.gesdisc.eosdis.nasa.gov:80/dods/GLDAS_NOAH025_M' %(username,password)
        url = 'https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_NOAH0125_M.002' #%(username,password)


        # Name the definition that will be used to obtain the data
        RetrieveData_fcn = RetrieveData_monthly

    # If none of the possible option are chosen
    else:
        raise KeyError("The input time interval is not supported")

    # Define IDs (latitude/longitude)
    yID = np.int16(np.array([np.ceil((latlim[0] - 25) * 8),
                             np.floor((latlim[1] - 25) * 8)]))
    xID = np.int16(np.array([np.floor((lonlim[0] + 125) * 8),
                             np.ceil((lonlim[1] + 125) * 8)]))

    # Check dates. If no dates are given, the max number of days is used.
    if not Startdate:
        Startdate = pd.Timestamp(sd_date)
    if not Enddate:
        Enddate = pd.Timestamp('Now')  # Should be much than available

    # Create all dates that will be calculated
    Dates = pd.date_range(Startdate, Enddate, freq=TimeFreq)

    # Create Waitbar
    if Waitbar == 1:
        import watertools.Functions.Random.WaitbarConsole as WaitbarConsole
        total_amount = len(Dates)
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    # Define the variable string name
    VarStr = VarInfo.names[Var]

    # Create one parameter with all the required arguments
    args = [path, url, Var, VarStr, VarInfo,
            TimeCase, xID, yID, lonlim, latlim, CaseParameters, username, password]

    # Pass variables to parallel function and run
    if not cores:
        for Date in Dates:
            RetrieveData_fcn(Date, args)
            if Waitbar == 1:
                amount += 1
                WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)
        results = True
    else:
        results = Parallel(n_jobs=cores)(delayed(RetrieveData_fcn)(Date, args)
                                         for Date in Dates)
    return results


def RetrieveData_hourly(Date, args):
    """
    This function retrieves NLDAS hourly data for a given date.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """

	# Open all the parameters
    [path, url, Var, VarStr, VarInfo, TimeCase, xID, yID, lonlim, latlim, CaseParameters, username, password] = args

	# Open variable info parameters
    VarFactor = VarInfo.factors[Var]


	# Loop over the periods
    for period in CaseParameters:

        hour = (int(period)-1)
        
        # Check whether the file already exist or the worldfile is
        # downloaded
        BasinDir = path + '/' + VarStr + '_NLDAS-NOAH_' + \
            VarInfo.units[Var] + '_hourly_' + Date.strftime('%Y.%m.%d') + \
            '_%02d00.tif' %(hour)

        if not os.path.isfile(BasinDir):

            # Reset the begin parameters for downloading
            downloaded = 0
            N=0

            while downloaded == 0:
                try:

                    # Define time
                    zID = int(((Date - pd.Timestamp("1979-1-2")).days) * 24) + (period - 1) - 1

                    # total URL
                    url_NLDAS = url + '.ascii?%s[%s][%s:1:%s][%s:1:%s]' %(Var,zID,yID[0],yID[1],xID[0],xID[1])

                    # open URL
                    try:
                        dataset = requests.get(url_NLDAS, allow_redirects=False,stream = True)
                    except:
                        from requests.packages.urllib3.exceptions import InsecureRequestWarning
                        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                        dataset = requests.get(url_NLDAS, allow_redirects=False,stream = True, verify = False)

                    try:
                        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True)
                    except:
                        from requests.packages.urllib3.exceptions import InsecureRequestWarning
                        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True, verify = False)

                    # download data (first save as text file)
                    pathtext = os.path.join(path, 'temp%s.txt' % zID)
                    z = open(pathtext, 'w')
                    z.write(get_dataset.text)
                    z.close()

                    # Open text file and remove header and footer
                    data_start = np.genfromtxt(pathtext,dtype = float,skip_header = 1,skip_footer = 6,delimiter=',')
                    data = data_start[:,1:]

                    # Add the VarFactor
                    if VarFactor < 0:
                        data = data + VarFactor
                    else:
                        data = data * VarFactor

                    # Set Nan value for values lower than -9999
                    data[data < -9999] = -9999

                    # Say that download was succesfull
                    downloaded = 1

                # If download was not succesfull
                except:
                    data=[]

                    # Try another time
                    N = N + 1

                    # Stop trying after 10 times
                    if N == 10:
                        print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                        downloaded = 1

            # define geo
            lonlimNLDAS = xID[0] * 0.125 - 125.0
            latlimNLDAS = (yID[1] + 1) * 0.125 + 25.0

            # Save to geotiff file
            geo = [lonlimNLDAS,0.125,0,latlimNLDAS,0,-0.125]
            DC.Save_as_tiff(name=BasinDir, data=np.flipud(data[:,:]), geo=geo, projection="WGS84")

            # Delete data and text file
            del data
            os.remove(pathtext)

    return True


def RetrieveData_daily(Date, args):
    """
    This function retrieves GLDAS daily data for a given date.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """

    # Open all the parameters
    [path, url, Var, VarStr, VarInfo,
     TimeCase, xID, yID, lonlim, latlim, CaseParameters, username, password] = args
    [selected, types] = CaseParameters

    # Reset the begin parameters for downloading
    downloaded = 0
    N = 0
    data_end = []

    # Open all variable info
    for T in types:
        if T == 'mean':
            VarStr = VarInfo.names[Var]
        else:
            VarStr = VarInfo.names[Var] + '-' + T

        # Check whether the file already exist or
        # the worldfile is downloaded
        BasinDir = path[T] + '/' + VarStr + '_NLDAS-NOAH_' + \
            VarInfo.units[Var] + '_daily_' + Date.strftime('%Y.%m.%d') + \
            '.tif'

        # Check if the outputfile already excists
        if not os.path.isfile(BasinDir):

            # Create the time dimension
            zID_start = int(((Date - pd.Timestamp("1979-1-2")).days) * 24) - 1
            zID_end = zID_start + 23


            # define total url
            url_NLDAS = url + '.ascii?%s[%s:1:%s][%s:1:%s][%s:1:%s]' %(Var,zID_start,zID_end,yID[0],yID[1],xID[0],xID[1])

            # if not downloaded try to download file
            while downloaded == 0:
                try:

                    # open URL
                    try:
                        dataset = requests.get(url_NLDAS, allow_redirects=False,stream = True)
                    except:
                        from requests.packages.urllib3.exceptions import InsecureRequestWarning
                        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                        dataset = requests.get(url_NLDAS, allow_redirects=False,stream = True, verify = False)
                    try:
                        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True)
                    except:
                        from requests.packages.urllib3.exceptions import InsecureRequestWarning
                        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True, verify = False)

                    # download data (first save as text file)
                    pathtext = os.path.join(path[T],'temp%s.txt' %str(zID_start))
                    z = open(pathtext,'w')
                    z.write(get_dataset.text)
                    z.close()

                    # Reshape data
                    datashape = [8,yID[1] - yID[0] + 1,xID[1] - xID[0] + 1]
                    data_start = np.genfromtxt(pathtext,dtype = float,skip_header = 1,skip_footer = 6,delimiter = ',')
                    data_list = np.asarray(data_start[:,1:])
                    data_end = np.resize(data_list,(8, datashape[1], datashape[2]))
                    os.remove(pathtext)

                    # Add the VarFactor
                    if VarInfo.factors[Var] < 0:
                        data_end[data_end != -9999] = data_end[data_end != -9999] + VarInfo.factors[Var]
                    else:
                        data_end[data_end != -9999] = data_end[data_end != -9999] * VarInfo.factors[Var]
                    if VarInfo.types[Var] == 'flux':
                        data_end = data_end * 24    
                    data_end[data_end < -9999] = -9999

                    # define geo
                    lonlimNLDAS = xID[0] * 0.125 - 125.0
                    latlimNLDAS = (yID[1] + 1) * 0.125 + 25.0
                
                    # Download was succesfull
                    downloaded = 1

                # If download was not succesfull
                except:

                    # Try another time
                    N = N + 1

                    # Stop trying after 10 times
                    if N == 10:
                        print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                        downloaded = 1

            try:
                # Save to geotiff file

                if T == 'mean':
                    data = np.flipud(np.mean(data_end, axis=0))
                if T == 'max':
                    data = np.flipud(np.max(data_end, axis=0))
                if T == 'min':
                    data = np.flipud(np.min(data_end, axis=0))

                geo = [lonlimNLDAS,0.125,0,latlimNLDAS,0,-0.125]
                DC.Save_as_tiff(name=BasinDir, data=data, geo=geo, projection="WGS84")

            except:
                print('NLDAS map from '+ Date.strftime('%Y-%m-%d') + ' is not created')

    return True


def RetrieveData_monthly(Date, args):
    """
    This function retrieves GLDAS monthly data for a given date.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """
    # Argument
    [path, url, Var, VarStr, VarInfo,
            TimeCase, xID, yID, lonlim, latlim, CaseParameters, username, password] = args

	# Open variable info parameters
    VarFactor = VarInfo.factors[Var]

    # Check whether the file already exist or the worldfile is downloaded
    BasinDir = path + '/' + VarStr + '_NLDAS-NOAH_' + \
        VarInfo.units[Var] + '_monthly_' + Date.strftime('%Y.%m.%d') + \
        '.tif'

    # Define month and year of current month
    Y = Date.year
    M = Date.month
    Mday = calendar.monthrange(Y, M)[1]

    # Check if the outputfile already excists
    if not os.path.isfile(BasinDir):

        # Reset the begin parameters for downloading
        downloaded = 0
        N=0

        # Create the time dimension
        zID = (Y - 1979) * 12 + (M - 1)

        # define total url
        url_NLDAS = url + '.ascii?%s[%s][%s:1:%s][%s:1:%s]' %(Var,zID,yID[0],yID[1],xID[0],xID[1])

        # if not downloaded try to download file
        while downloaded == 0:
            try:

                # open URL
                try:
                    dataset = requests.get(url_NLDAS, allow_redirects=False,stream = True)
                except:
                    from requests.packages.urllib3.exceptions import InsecureRequestWarning
                    requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                    dataset = requests.get(url_NLDAS, allow_redirects=False,stream = True, verify = False)
                try:
                    get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True)
                except:
                    from requests.packages.urllib3.exceptions import InsecureRequestWarning
                    requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                    get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True, verify = False)

                # download data (first save as text file)
                pathtext = os.path.join(path,'temp%s.txt' %str(zID))
                z = open(pathtext,'w')
                z.write(get_dataset.text)
                z.close()

                # Open text file and remove header and footer
                data_start = np.genfromtxt(pathtext,dtype = float,skip_header = 1,skip_footer = 6,delimiter=',')
                data = data_start[:,1:]

                # Add the VarFactor
                if VarFactor < 0:
                    data = data + VarFactor
                else:
                    data = data * VarFactor
                if VarInfo.types[Var] == 'flux':
                    data = data * Mday  
                    
                # Set Nan value for values lower than -9999
                data[data < -9999] = -9999

                # Say that download was succesfull
                downloaded = 1

            # If download was not succesfull
            except:
                data=[]

                # Try another time
                N = N + 1

                # Stop trying after 10 times
                if N == 10:
                    print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                    downloaded = 1

            # define geo
            lonlimNLDAS = xID[0] * 0.125 - 125.0
            latlimNLDAS = (yID[1] + 1) * 0.125 + 25.0

            # Save to geotiff file
            geo = [lonlimNLDAS,0.125,0,latlimNLDAS,0,-0.125]
            DC.Save_as_tiff(name=BasinDir, data=np.flipud(data[:,:]), geo=geo, projection="WGS84")


            # Delete data and text file
            del data
            #os.remove(pathtext)

    return True

class VariablesInfo:
    """
    This class contains the information about the GLDAS variables
    """
    names = {'acondsfc': 'Aero_Conductance',
             'albdosfc': 'Albedo',
             'arainsfc': 'P',
             'asnowsfc': 'Snow',
             'avsftsfc': 'Surface_Temp',
             'bgrunsfc': 'Subsurface_Runoff',
             'ccondsfc': 'Canopy_Conductance',
             'cnwatsfc': 'Plant_Canopy_Sur_Water',
             'dlwrfsfc': 'DLWR',
             'dswrfsfc': 'DSWR',
             'evbssfc': 'Direct_Evap_Bare_Soil',             
             'evcwsfc': 'Canopy_Evaporation',
             'evpsfc': 'ETtot',
             'gfluxsfc': 'Ground_Heat_Flux',
             'laisfc': 'LAI',
             'lhtflsfc': 'LE',
             'lsoil0_10cm': 'Liquid_Soil_Moisture_Content_0_10cm',             
             'lsoil10_40cm': 'Liquid_Soil_Moisture_Content_10_40cm',
             'lsoil40_100cm': 'Liquid_Soil_Moisture_Content_40_100cm',
             'lsoil100_200cm': 'Liquid_Soil_Moisture_Content_100_200cm',
             'mstav0_100cm': 'Moisture_Availability_0_100cm',
             'mstav0_200cm': 'Moisture_Availability_0_200cm',
             'nlwrssfc': 'NLWR',             
             'nswrssfc': 'NSWR',
             'pevprsfc': 'ETpot',
             'rcqsfc': 'Hum_in_Canopy_Conductance',
             'rcssfc': 'Solar_in_Canopy_Conductance',
             'rcsolsfc': 'Soil_Moisture_in_Canopy_Conductance',
             'rctsfc': 'Temp_in_Canopy_Conductance',             
             'rsmacrsfc': 'Rel_Soil_Moisture_Availabitily',
             'rsminsfc': 'Min_Stomatal_Resistance',
             'rzsmgnd': 'Root_Zone_Soil_Moisture',
             'sbsnosfc': 'Sublimation',
             'shtflsfc': 'S',
             'snodsfc': 'Snow_Depth',   
             'snohfsfc':'Snow_Phase_Heat_Flux',
             'snomsfc': 'Snow_Melt',
             'snowcsfc': 'Snow_Cover',
             'soilm0_10cm': 'Soil_Moisture_Content_L1',
             'soilm0_100cm': 'Soil_Moisture_Content_Top',
             'soilm0_200cm': 'Soil_Moisture_Content_Total',
             'soilm10_40cm': 'Soil_Moisture_Content_L2',             
             'soilm40_100cm': 'Soil_Moisture_Content_L3',
             'soilm100_200cm': 'Soil_Moisture_Content_L4',
             'ssrunsfc': 'Surface_Runoff',
             'transsfc': 'T',
             'tsoil0_10cm': 'Soil_Temp_L1',
             'tsoil10_40cm': 'Soil_Temp_L2',
             'tsoil40_100cm': 'Soil_Temp_L3',
             'tsoil100_200cm': 'Soil_Temp_L4',
             'vegsfc': 'Veg',
             'weasdsfc': 'Acc_Snow_Water_Equivalent'}
    
    descriptions ={'acondsfc': 'aerodynamic conductance [m/s] ',
             'albdosfc': 'albedo [%] ',
             'arainsfc': 'rainfall (unfrozen precipitation) [kg/m^2] ',
             'asnowsfc': 'snowfall (frozen precipitation) [kg/m^2] ',
             'avsftsfc': 'average surface skin temperature [k] ',
             'bgrunsfc': 'subsurface runoff (baseflow) [kg/m^2] ',
             'ccondsfc': 'canopy conductance [m/s] ',
             'cnwatsfc': 'plant canopy surface water [kg/m^2] ',
             'dlwrfsfc': 'longwave radiation flux downwards (surface) [w/m^2] ',
             'dswrfsfc': 'shortwave radiation flux downwards (surface) [w/m^2] ',
             'evbssfc': 'direct evaporation from bare soil [w/m^2] ',             
             'evcwsfc': 'canopy water evaporation [w/m^2] ',
             'evpsfc': 'total evapotranspiration [kg/m^2] ',
             'gfluxsfc': 'ground heat flux [w/m^2] ',
             'laisfc': 'leaf area index (0-9) [unitless] ',
             'lhtflsfc': 'latent heat flux [w/m^2] ',
             'lsoil0_10cm': '0-10 cm liquid soil moisture content (non-frozen) [kg/m^2] ',             
             'lsoil10_40cm': '10-40 cm liquid soil moisture content (non-frozen) [kg/m^2] ',
             'lsoil40_100cm': '40-100 cm liquid soil moisture content (non-frozen) [kg/m^2] ',
             'lsoil100_200cm': '100-200 cm liquid soil moisture content (non-frozen) [kg/m^2] ',
             'mstav0_100cm': '0-100 cm moisture availability [%] ',
             'mstav0_200cm': ' 0-200 cm moisture availability [%] ',
             'nlwrssfc': 'longwave radiation flux net (surface) [w/m^2] ',             
             'nswrssfc': 'shortwave radiation flux net (surface) [w/m^2] ',
             'pevprsfc': 'potential evaporation rate [w/m^2] ',
             'rcqsfc': 'humidity parameter in canopy conductance [fraction] ',
             'rcssfc': 'solar parameter in canopy conductance [fraction] ',
             'rcsolsfc': 'soil moisture parameter in canopy conductance [fraction] ',
             'rctsfc': 'temperature parameter in canopy conductance [fraction] ',             
             'rsmacrsfc': 'relative soil moisture availability control factor [0-1] ',
             'rsminsfc': 'minimal stomatal resistance [s/m] ',
             'rzsmgnd': 'root zone soil moisture [kg/m^2] ',
             'sbsnosfc': 'sublimation (evaporation from snow) [w/m^2] ',
             'shtflsfc': 'sensible heat flux [w/m^2] ',
             'snodsfc': 'snow depth [m] ',
             'snohfsfc': 'snow phase-change heat flux [w/m^2] ',
             'snomsfc': ' snow melt [kg/m^2]  ',
             'snowcsfc': 'snow cover [fraction]  ',
             'soilm0_10cm': '0-10 cm layer 1 soil moisture content [kg/m^2] ',
             'soilm0_100cm': '0-100 cm top 1 meter soil moisture content [kg/m^2] ',
             'soilm0_200cm': '0-200 cm total column soil moisture content [kg/m^2] ',
             'soilm10_40cm': '10-40 cm layer 2 soil moisture content [kg/m^2] ',             
             'soilm40_100cm': '40-100 cm layer 3 soil moisture content [kg/m^2] ',
             'soilm100_200cm': '100-200 cm layer 4 soil moisture content [kg/m^2]',
             'ssrunsfc': 'surface runoff (non-infiltrating) [kg/m^2] ',
             'transsfc': 'transpiration [w/m^2] ',
             'tsoil0_10cm': '0-10 cm soil temperature [k] ',
             'tsoil10_40cm': '10-40 cm soil temperature [k] ',
             'tsoil40_100cm': '40-100 cm soil temperature [k] ',
             'tsoil100_200cm': '100-200 cm soil temperature [k] ',
             'vegsfc': 'vegetation [fraction] ',
             'weasdsfc': 'accumulated snow water-equivalent [kg/m^2] '}
    factors = {'acondsfc': 1,
             'albdosfc': 1,
             'arainsfc': 1,
             'asnowsfc': 1,
             'avsftsfc': -273.15,
             'bgrunsfc': 1,
             'ccondsfc': 1,
             'cnwatsfc': 1,
             'dlwrfsfc': 1,
             'dswrfsfc': 1,
             'evbssfc': 1,           
             'evcwsfc': 1,
             'evpsfc': 1,
             'gfluxsfc': 1,
             'laisfc': 1,
             'lhtflsfc': 1,
             'lsoil0_10cm': 1,          
             'lsoil10_40cm': 1,
             'lsoil40_100cm': 1,
             'lsoil100_200cm': 1,
             'mstav0_100cm': 1,
             'mstav0_200cm': 1,
             'nlwrssfc': 1,           
             'nswrssfc': 1,
             'pevprsfc': 1,
             'rcqsfc': 1,
             'rcssfc': 1,
             'rcsolsfc': 1,
             'rctsfc': 1,         
             'rsmacrsfc': 1,
             'rsminsfc': 1,
             'rzsmgnd': 1,
             'sbsnosfc': 1,
             'shtflsfc': 1,
             'snodsfc': 1,
             'snohfsfc': 1,
             'snomsfc': 1,
             'snowcsfc': 1,
             'soilm0_10cm': 1,
             'soilm0_100cm': 1,
             'soilm0_200cm': 1,
             'soilm10_40cm': 1,           
             'soilm40_100cm': 1,
             'soilm100_200cm': 1,
             'ssrunsfc': 1,
             'transsfc': 1,
             'tsoil0_10cm': -273.15,
             'tsoil10_40cm': -273.15,
             'tsoil40_100cm': -273.15,
             'tsoil100_200cm': -273.15,
             'vegsfc': 1,
             'weasdsfc': 1}
    
    types = {'acondsfc': 'state',
             'albdosfc': 'state',
             'arainsfc': 'flux',
             'asnowsfc': 'flux',
             'avsftsfc': 'state',
             'bgrunsfc': 'flux',
             'ccondsfc': 'state',
             'cnwatsfc': 'flux',
             'dlwrfsfc': 'state',
             'dswrfsfc': 'state',
             'evbssfc': 'flux',            
             'evcwsfc': 'flux',
             'evpsfc': 'flux',
             'gfluxsfc': 'flux',
             'laisfc': 'state',
             'lhtflsfc': 'flux',
             'lsoil0_10cm': 'state',             
             'lsoil10_40cm': 'state',
             'lsoil40_100cm': 'state',
             'lsoil100_200cm': 'state',
             'mstav0_100cm': 'state',
             'mstav0_200cm': 'state',
             'nlwrssfc': 'state',             
             'nswrssfc': 'state',
             'pevprsfc': 'flux',
             'rcqsfc': 'state',
             'rcssfc': 'state',
             'rcsolsfc': 'state',
             'rctsfc': 'state',             
             'rsmacrsfc': 'state',
             'rsminsfc': 'state',
             'rzsmgnd': 'state',
             'sbsnosfc': 'flux',
             'shtflsfc': 'flux',
             'snodsfc': 'state',
             'snohfsfc': 'flux',
             'snomsfc': 'flux',
             'snowcsfc': 'state',
             'soilm0_10cm': 'state',
             'soilm0_100cm': 'state',
             'soilm0_200cm': 'state',
             'soilm10_40cm': 'state',             
             'soilm40_100cm': 'state',
             'soilm100_200cm': 'state',
             'ssrunsfc': 'flux',
             'transsfc': 'flux',
             'tsoil0_10cm': 'state',
             'tsoil10_40cm': 'state',
             'tsoil40_100cm': 'state',
             'tsoil100_200cm': 'state',
             'vegsfc': 'state',
             'weasdsfc': 'flux'}

    def __init__(self, step):
        if step == 'hourly':
            self.units = {'acondsfc': 'm-s',
             'albdosfc': 'Percentage',
             'arainsfc': 'mm-hour-1',
             'asnowsfc': 'mm-hour-1',
             'avsftsfc': 'C',
             'bgrunsfc': 'mm-hour-1',
             'ccondsfc': 'm-s-1',
             'cnwatsfc': 'mm-hour-1',
             'dlwrfsfc': 'w-m-2',
             'dswrfsfc': 'w-m-2',
             'evbssfc': 'w-m-2',             
             'evcwsfc': 'w-m-2',
             'evpsfc': 'mm-hour-1',
             'gfluxsfc': 'w-m-2',
             'laisfc': '-',
             'lhtflsfc': 'w-m-2',
             'lsoil0_10cm': 'mm-hour-1',             
             'lsoil10_40cm': 'mm-hour-1',
             'lsoil40_100cm': 'mm-hour-1',
             'lsoil100_200cm': 'mm-hour-1',
             'mstav0_100cm': 'Percentage',
             'mstav0_200cm': 'Percentage',
             'nlwrssfc': 'w-m-2',             
             'nswrssfc': 'w-m-2',
             'pevprsfc': 'w-m-2',
             'rcqsfc': '-',
             'rcssfc': '-',
             'rcsolsfc': '-',
             'rctsfc': '-',             
             'rsmacrsfc': '-',
             'rsminsfc': 's-m-1',
             'rzsmgnd': 'mm-hour-1',
             'sbsnosfc': 'w-m-2',
             'shtflsfc': 'w-m-2',
             'snodsfc': 'm',
             'snohfsfc': 'w-m-2',
             'snomsfc': 'mm-hour-1',
             'snowcsfc': '-',
             'soilm0_10cm': 'mm',
             'soilm0_100cm': 'mm',
             'soilm0_200cm': 'mm',
             'soilm10_40cm': 'mm',             
             'soilm40_100cm': 'mm',
             'soilm100_200cm': 'm',
             'ssrunsfc': 'mm-hour-1',
             'transsfc': 'w-m-2',
             'tsoil0_10cm': 'C',
             'tsoil10_40cm': 'C',
             'tsoil40_100cm': 'C',
             'tsoil100_200cm': 'C',
             'vegsfc': '-',
             'weasdsfc': 'mm'}
            
            
        elif step == 'daily':
            self.units = {'acondsfc': 'm-s',
             'albdosfc': 'Percentage',
             'arainsfc': 'mm-day-1',
             'asnowsfc': 'mm-day-1',
             'avsftsfc': 'C',
             'bgrunsfc': 'mm-day-1',
             'ccondsfc': 'm-s-1',
             'cnwatsfc': 'mm-day-1',
             'dlwrfsfc': 'w-m-2',
             'dswrfsfc': 'w-m-2',
             'evbssfc': 'w-m-2',             
             'evcwsfc': 'w-m-2',
             'evpsfc': 'mm-day-1',
             'gfluxsfc': 'w-m-2',
             'laisfc': '-',
             'lhtflsfc': 'w-m-2',
             'lsoil0_10cm': 'mm-day-1',             
             'lsoil10_40cm': 'mm-day-1',
             'lsoil40_100cm': 'mm-day-1',
             'lsoil100_200cm': 'mm-day-1',
             'mstav0_100cm': 'Percentage',
             'mstav0_200cm': 'Percentage',
             'nlwrssfc': 'w-m-2',             
             'nswrssfc': 'w-m-2',
             'pevprsfc': 'w-m-2',
             'rcqsfc': '-',
             'rcssfc': '-',
             'rcsolsfc': '-',
             'rctsfc': '-',             
             'rsmacrsfc': '-',
             'rsminsfc': 's-m-1',
             'rzsmgnd': 'mm-day-1',
             'sbsnosfc': 'w-m-2',
             'shtflsfc': 'w-m-2',
             'snodsfc': 'm',
             'snohfsfc': 'w-m-2',
             'snomsfc': 'mm-day-1',
             'snowcsfc': '-',
             'soilm0_10cm': 'mm',
             'soilm0_100cm': 'mm',
             'soilm0_200cm': 'mm',
             'soilm10_40cm': 'mm',             
             'soilm40_100cm': 'mm',
             'soilm100_200cm': 'm',
             'ssrunsfc': 'mm-day-1',
             'transsfc': 'w-m-2',
             'tsoil0_10cm': 'C',
             'tsoil10_40cm': 'C',
             'tsoil40_100cm': 'C',
             'tsoil100_200cm': 'C',
             'vegsfc': '-',
             'weasdsfc': 'mm'}
        elif step == 'monthly':
            self.units = {'acondsfc': 'm-s',
             'albdosfc': 'Percentage',
             'arainsfc': 'mm-month-1',
             'asnowsfc': 'mm-month-1',
             'avsftsfc': 'C',
             'bgrunsfc': 'mm-month-1',
             'ccondsfc': 'm-s-1',
             'cnwatsfc': 'mm-month-1',
             'dlwrfsfc': 'w-m-2',
             'dswrfsfc': 'w-m-2',
             'evbssfc': 'w-m-2',             
             'evcwsfc': 'w-m-2',
             'evpsfc': 'mm-month-1',
             'gfluxsfc': 'w-m-2',
             'laisfc': '-',
             'lhtflsfc': 'w-m-2',
             'lsoil0_10cm': 'mm-month-1',             
             'lsoil10_40cm': 'mm-month-1',
             'lsoil40_100cm': 'mm-month-1',
             'lsoil100_200cm': 'mm-month-1',
             'mstav0_100cm': 'Percentage',
             'mstav0_200cm': 'Percentage',
             'nlwrssfc': 'w-m-2',             
             'nswrssfc': 'w-m-2',
             'pevprsfc': 'w-m-2',
             'rcqsfc': '-',
             'rcssfc': '-',
             'rcsolsfc': '-',
             'rctsfc': '-',             
             'rsmacrsfc': '-',
             'rsminsfc': 's-m-1',
             'rzsmgnd': 'mm-month-1',
             'sbsnosfc': 'w-m-2',
             'shtflsfc': 'w-m-2',
             'snodsfc': 'm',
             'snohfsfc': 'w-m-2',
             'snomsfc': 'mm-month-1',
             'snowcsfc': '-',
             'soilm0_10cm': 'mm',
             'soilm0_100cm': 'mm',
             'soilm0_200cm': 'mm',
             'soilm10_40cm': 'mm',             
             'soilm40_100cm': 'mm',
             'soilm100_200cm': 'm',
             'ssrunsfc': 'mm-month-1',
             'transsfc': 'w-m-2',
             'tsoil0_10cm': 'C',
             'tsoil10_40cm': 'C',
             'tsoil40_100cm': 'C',
             'tsoil100_200cm': 'C',
             'vegsfc': '-',
             'weasdsfc': 'mm'}
        else:
            raise KeyError("The input time step is not supported")
