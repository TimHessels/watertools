# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:27:36 2016

@author: tih
"""


def Accounts(Type=None):

    User_Pass = {
     'NASA': ['TimHessels', 'WAteam1!'],
     'GLEAM': ['gleamuser', 'hi_GLEAMv32#HiL?'],
     'FTP_WA': ['THessels', 'painole_2016!'],
     'MSWEP': ['THessels', '9rvkvkG6z'],
     'Copernicus': ['timhessels91', 'EelsBloc91!2'],  #https://land.copernicus.vgt.vito.be/PDF/
     'VITO': ['TimHessels91WA', 'EelsBloc91!2']}     #https://www.vito-eodata.be/PDF/datapool/
	 
    Selected_Path = User_Pass[Type]

    return(Selected_Path)
