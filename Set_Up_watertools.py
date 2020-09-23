# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 08:28:08 2020

@author: timhe
"""
import os
import json

from cryptography.fernet import Fernet

def create_key():
    
    key_filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), "wa_key.txt")
    
    if os.path.exists(key_filename):
        os.remove(key_filename)
    
    f = open(key_filename,"w+")
    f.write(str(Fernet.generate_key().decode("utf-8")))
    f.close()

    return()

def set_up_account(account="ALL"):
    
    key_filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), "wa_key.txt")
    json_filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), "keys.json")
    
    if not os.path.exists(key_filename):
        create_key()
    
    f = open(key_filename,"r")
    key = f.read()
    f.close()
    
    cipher_suite = Fernet(key.encode('utf-8'))

    if os.path.exists(json_filename):
    
        # open json file  
        with open(json_filename) as f:
            datastore = f.read()
        obj = json.loads(datastore)      
        f.close()
    else:
        obj = {}

    if account == "ALL":
        Servers = ["NASA", "GLEAM", "FTP_WA", "MSWEP", "Copernicus", "VITO", "WAPOR"]
    else:
        Servers = [account]
    
    for Server in Servers:
        
        if Server != "WAPOR":
            account_name = input("Type in your account username for %s" %Server)
            pwd = input("Type in your password for %s" %Server)        
            
            account_name_crypt = cipher_suite.encrypt(("%s" %account_name).encode('utf-8'))
            pwd_crypt = cipher_suite.encrypt(("%s" %pwd).encode('utf-8'))
            obj[Server] = ([str(account_name_crypt.decode("utf-8")), str(pwd_crypt.decode("utf-8"))])
        if  Server == "WAPOR":
            API_Key = input("Type in your API key for %s" %Server)
            API_Key_crypt = cipher_suite.encrypt(("%s" %API_Key).encode('utf-8'))
            obj[Server] = [str(API_Key_crypt.decode("utf-8"))]
            
    # save extent in task
    with open(json_filename, 'w') as outfile:
        json.dump(obj, outfile)        
    
    return()
    
