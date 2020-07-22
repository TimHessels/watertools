# -*- coding: utf-8 -*-2019
"""
Created on Thu Mar 17 09:27:36 2016

@author: tih
"""


def Accounts(Type=None):

    User_Pass = {
     'NASA': ['TimHessels', 'WAteam1!'],
     'GLEAM': ['gleamuser', 'v33_GLEAM2019#aw'],
     'FTP_WA': ['THessels', 'painole_2016!'],
     'MSWEP': ['THessels', '9rvkvkG6z'],
     'Copernicus': ['timhessels91', 'EelsBloc91!2'],  #https://land.copernicus.vgt.vito.be/PDF/
     'VITO': ['TimHessels91WA', 'EelsBloc91!2'],#https://www.vito-eodata.be/PDF/datapool/
     'WAPOR': ['02d7f2b8e00f26f2212453ade00718ccc9258c4c4166ef82d4eec2aab58923ccc12aeaaac4f150b8']}     
	 
    Selected_Path = User_Pass[Type]

    return(Selected_Path)
