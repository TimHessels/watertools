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
        path = os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS_FORA',
                            TimeCase, Var)

        if not os.path.exists(path):
            os.makedirs(path)

        # Startdate if not defined
        sd_date = '1979-01-02'

        # Define Time frequency
        TimeFreq = 'D'

        # Define URL by using personal account
        #url = 'http://%s:%s@hydro1.gesdisc.eosdis.nasa.gov:80/dods/GLDAS_NOAH025SUBP_3H' %(username,password)
        url = 'https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_FORA0125_H.002' #%(username,password)

        # Name the definition that will be used to obtain the data
        RetrieveData_fcn = RetrieveData_hourly

    # Set required data for the daily option
    elif TimeCase == 'daily':

        # seperate the daily case parameters
        SumMean, Min, Max = CaseParameters

        # Define output folder and create this one if not exists
        path = {'mean': os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS_FORA',
                                     TimeCase, Var, 'mean'),
                'min': os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS_FORA',
                                    TimeCase, Var, 'min'),
                'max': os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS_FORA',
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
        url = 'https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_FORA0125_H.002' #%(username,password)

        # Name the definition that will be used to obtain the data
        RetrieveData_fcn = RetrieveData_daily

    # Set required data for the monthly option
    elif TimeCase == 'monthly':

        # Define output folder and create this one if not exists
        path = os.path.join(Dir, 'Weather_Data', 'Model', 'NLDAS_FORA',
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
        url = 'https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_FORA0125_M.002' #%(username,password)


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
        BasinDir = path + '/' + VarStr + '_NLDAS-FORA_' + \
            VarInfo.units[Var] + '_hourly_' + Date.strftime('%Y.%m.%d') + \
            '_%02d00.tif' %(hour)
        
        if not os.path.isfile(BasinDir):

            # Reset the begin parameters for downloading
            downloaded = 0
            N=0

            while downloaded == 0:
                try:

                    # Define time
                    zID = int(((Date - pd.Timestamp("1979-1-2")).days) * 24) + (period - 1) - 13

                    # total URL
                    url_NLDAS = url + '.ascii?%s[%s][%s:1:%s][%s:1:%s]' %(Var,zID,yID[0],yID[1],xID[0],xID[1])
                    #print(url_NLDAS)
                    # open URL
                    try:
                        dataset = requests.get(url_NLDAS, allow_redirects=False,stream = True)
                    except:
                        from requests.packages.urllib3.exceptions import InsecureRequestWarning
                        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                        dataset = requests.get(url_NLDAS, allow_redirects=False,stream = True, verify = False)
                    
                    '''
                    try:
                        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True)
                    except:
                        from requests.packages.urllib3.exceptions import InsecureRequestWarning
                        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True, verify = False)
                    '''
                    # download data (first save as text file)
                    pathtext = os.path.join(path, 'temp%s.txt' % zID)

                    z = open(pathtext, 'w')
                    z.write(dataset.text)
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
                    if N == 5:
                        print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                        print(url_NLDAS)
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
        BasinDir = path[T] + '/' + VarStr + '_NLDAS-FORA_' + \
            VarInfo.units[Var] + '_daily_' + Date.strftime('%Y.%m.%d') + \
            '.tif'

        # Check if the outputfile already excists
        if not os.path.isfile(BasinDir):

            # Create the time dimension
            zID_start = int(((Date - pd.Timestamp("1979-1-2")).days) * 24) - 13
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
                        
                    '''
                    try:
                        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True)
                    except:
                        from requests.packages.urllib3.exceptions import InsecureRequestWarning
                        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                        get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True, verify = False)
                    '''
                    
                    # download data (first save as text file)
                    pathtext = os.path.join(path[T],'temp%s.txt' %str(zID_start))
                    z = open(pathtext,'w')
                    z.write(dataset.text)
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
                    if N == 5:
                        print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                        print(url_NLDAS)
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
    BasinDir = path + '/' + VarStr + '_NLDAS-FORA_' + \
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
                    
                '''
                try:
                    get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True)
                except:
                    from requests.packages.urllib3.exceptions import InsecureRequestWarning
                    requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
                    get_dataset = requests.get(dataset.headers['location'], auth = (username,password),stream = True, verify = False)
                '''
                
                # download data (first save as text file)
                pathtext = os.path.join(path,'temp%s.txt' %str(zID))
                z = open(pathtext,'w')
                z.write(dataset.text)
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
                if N == 5:
                    print('Data from ' + Date.strftime('%Y-%m-%d') + ' is not available')
                    print(url_NLDAS)
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
    names = {'apcpsfc': 'P',
             'convfracsfc': 'FracConvP',
             'dlwrfsfc': 'LWdown',
             'dswrfsfc': 'SWdown',
             'pevapsfc': 'ETpot',
             'pressfc': 'P',
             'spfh2m': 'Hum',
             'tmp2m': 'T2m',
             'ugrd10m': 'Wind10mZonal',
             'vgrd10m': 'Wind10mMeridional'
             }
    descriptions = {'apcpsfc': 'precipitation hourly total [kg/m^2]',
                    'convfracsfc': 'fraction of total precipitation that is convective [unitless]',
                    'dlwrfsfc': ' lw radiation flux downwards (surface) [w/m^2]',
                    'dswrfsfc': 'sw radiation flux downwards (surface) [w/m^2]',
                    'pevapsfc': 'potential evaporation [kg/m^2] ',
                    'pressfc': 'surface pressure [pa] ',
                    'spfh2m': '2-m above ground specific humidity [kg/kg]',
                    'tmp2m': '2-m above ground temperature [k] ',
                    'ugrd10m': '10-m above ground zonal wind speed [m/s]',
                    'vgrd10m': '10-m above ground meridional wind speed [m/s]'
                    }
    factors = {'apcpsfc': 1,
               'convfracsfc': 1,
               'dlwrfsfc': 1,
               'dswrfsfc': 1,
               'pevapsfc': 1,
               'pressfc': 0.001,
               'spfh2m': 1,
               'tmp2m': -273.15,
               'ugrd10m': 1,
               'vgrd10m': 1
             }
    types = {'apcpsfc': 'flux',
               'convfracsfc': 'state',
               'dlwrfsfc': 'state',
               'dswrfsfc': 'state',
               'pevapsfc': 'flux',
               'pressfc': 'state',
               'spfh2m': 'state',
               'tmp2m': 'state',
               'ugrd10m': 'state',
               'vgrd10m': 'state'
             }

    def __init__(self, step):
        if step == 'hourly':
            self.units = {'apcpsfc': 'mm-hour',
                          'convfracsfc': 'fraction',
                          'dlwrfsfc': 'W-m-2',
                          'dswrfsfc': 'W-m-2',
                          'pevapsfc': 'mm-hour',
                          'pressfc': 'kpa',
                          'spfh2m': 'kg-kg',
                          'tmp2m': 'C',
                          'ugrd10m': 'm-s-1',
                          'vgrd10m': 'm-s-1'
                          }
            
            
        elif step == 'daily':
            self.units = {'apcpsfc': 'mm-day',
                          'convfracsfc': 'fraction',
                          'dlwrfsfc': 'W-m-2',
                          'dswrfsfc': 'W-m-2',
                          'pevapsfc': 'mm-day',
                          'pressfc': 'kpa',
                          'spfh2m': 'kg-kg',
                          'tmp2m': 'C',
                          'ugrd10m': 'm-s-1',
                          'vgrd10m': 'm-s-1'
                          }
        elif step == 'monthly':
            self.units = {'apcpsfc': 'mm-month',
                          'convfracsfc': 'fraction',
                          'dlwrfsfc': 'W-m-2',
                          'dswrfsfc': 'W-m-2',
                          'pevapsfc': 'mm-month',
                          'pressfc': 'kpa',
                          'spfh2m': 'kg-kg',
                          'tmp2m': 'C',
                          'ugrd10m': 'm-s-1',
                          'vgrd10m': 'm-s-1'
                          }
        else:
            raise KeyError("The input time step is not supported")
