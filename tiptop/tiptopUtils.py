import os
import itertools
from datetime import datetime
from configparser import ConfigParser
import yaml

import numpy as np
from scipy.interpolate import interp1d

from matplotlib import rc

from mastsel import gpuEnabled as gpuMastsel
from p3.aoSystem import gpuEnabled as gpuP3

def arrayP3toMastsel(v):
    if (gpuP3 and gpuMastsel) or (not gpuP3 and not gpuMastsel):
        return v
    elif not gpuP3 and gpuMastsel:
        return cp.asarray(v)
    elif gpuP3 and not gpuMastsel: 
        return v.get()
    
def cpuArray(v):
    if isinstance(v,np.ndarray):
        return v
    else:
        return v.get()

MAX_VALUE_CHARS = 80
APPEND_TOKEN = '&&&'

def add_hdr_keyword(hdr, key_primary, key_secondary, val, iii=None, jjj=None):
    val_string = str(val)
    key = 'HIERARCH '+ key_primary +' '+ key_secondary
    if not iii is None:
        key += ' '+str(iii)    
    if not jjj is None:
        key +=  ' '+str(jjj)    
    margin = 4
    key = key
    current_val_string = val_string
    if len(key) + margin > MAX_VALUE_CHARS:
        print("Error, keywork is not acceptable due to string length.")
        return
    while not len(key) + 1 + len(current_val_string) + margin < MAX_VALUE_CHARS:
        max_char_index = MAX_VALUE_CHARS-len(key)-1-len(APPEND_TOKEN)-margin
        hdr[key+'+'] = current_val_string[:max_char_index]+APPEND_TOKEN        
        current_val_string = current_val_string[max_char_index:]        
    hdr[key] = current_val_string

def hdr2map(hdr):
    hdr_keys = list(hdr.keys())
    my_data_map = {}
    curr_value = ''
    curr_key = ''
    for key in hdr_keys:
        separator_index = key.find(' ')
        if separator_index > 0:            
            section = key[0:separator_index]
            if not section in my_data_map:
                my_data_map[section] = {}
            ext_indx = key.find('+')
            curr_key += key[separator_index+1:ext_indx]
            val_str = str(hdr[key])
            val_last_index = val_str.find(APPEND_TOKEN)            
            curr_value += val_str[:val_last_index]            
            if ext_indx==-1:
                my_data_map[section].update({curr_key:curr_value})
                curr_value = ''
                curr_key = ''                   
    return my_data_map

