# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:31:37 2020

@author: timhe
"""
import os
import sys
import json
import watertools

from cryptography.fernet import Fernet

def GET(server):
    
    path = os.path.dirname(watertools.__file__)
    
    key_file = os.path.join(path, "wa_key.txt")
    
    if not os.path.exists(key_file):
        sys.exit("Run watertools.Set_Up_watertools.create_key()")   
                 
    json_file = os.path.join(path, "keys.json")
    
    if not os.path.exists(json_file):
        sys.exit("Run watertools.Set_Up_watertools.set_up_account()")
        
    f = open(key_file,"r")
    key = f.read()
    f.close()
    
    cipher_suite = Fernet(key.encode('utf-8'))    
    
    # open json file  
    with open(json_file) as f:
        datastore = f.read()
    obj = json.loads(datastore)      
    f.close()
    
    if not server == "WAPOR":
        username_crypt, pwd_crypt = obj[server]
        username = cipher_suite.decrypt(username_crypt.encode("utf-8"))
        pwd = cipher_suite.decrypt(pwd_crypt.encode("utf-8"))        
    else:
        username_crypt = obj[server]
        username = cipher_suite.decrypt(username_crypt[0].encode("utf-8"))
        pwd = b''
    
    return(str(username.decode("utf-8")), str(pwd.decode("utf-8")))