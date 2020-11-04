# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:14:43 2020

@author: timhe
"""
import numpy as np

def All_Arrays(y, x):
    
    return_dict = dict()
    paras = dir()
    for para in paras:
        if type(para) == str and para[0] != "_":
            print(para)
            A = eval(para)
            
            if type(A) == np.ndarray:
                if A.shape[0]>y and A.shape[1]>x:
                    return_dict[para] = A[y,x]
                    return_dict["Shape_%s"%para] = [A.shape]
                           
    return(return_dict)
                
                
        
        