import os
import itertools
from datetime import datetime
from configparser import ConfigParser
import yaml

import numpy as np
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib import ticker
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
    if isinstance(v,np.ndarray) or isinstance(v, np.float64):
        return v
    else:
        return v.get()

MAX_VALUE_CHARS = 80
APPEND_TOKEN = '&&&'

def add_hdr_keyword(hdr, key_primary, key_secondary, val, iii=None, jjj=None):
    '''
    This functions add an element of the parmaters dictionary into the fits file header
    '''
    val_string = str(val)
    key = 'HIERARCH '+ key_primary +' '+ key_secondary
    if iii != None:
        key += ' '+str(iii)    
    if jjj != None:
        key += ' '+str(jjj)    
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
    '''
    Conversion of a fits file header into a dictionary
    '''
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

def plot_directions(parser, ticks_interval=5, labels=None):
    '''
    Polar plot with science and GS (sources_HO and sources_LO) directions

    :param parser: required, parameters object
    :type parser: configparser object
    :param ticks_interval: optional default=5, size of interval for ticks in the figure 
    :type ticks_interval: int
    :param labels: optional default=None, list of strings to be plotted next to the science sources
    :type labels: list
    
    :return: fig, ax
    :rtype: objects
    '''
    # SCIENCE
    if 'PSF_DIRECTIONS' in parser.sections(): # For retro-compatibility
        th_sci = np.array(eval(parser['PSF_DIRECTIONS']['ScienceAzimuth']))
        rr_sci = np.array(eval(parser['PSF_DIRECTIONS']['ScienceZenith']))
    else:
        th_sci = np.array(eval(parser['sources_science']['Azimuth']))
        rr_sci = np.array(eval(parser['sources_science']['Zenith']))
    th_sci = th_sci/180*np.pi
    fig = plt.figure('SOURCES DIRECTIONS', figsize=(6,6))
    ax = fig.add_subplot(111, polar=True)
    ax.tick_params(labelsize=10)
    ax.set_rlabel_position(225) # theta position of the radius labels
    ax.scatter(th_sci, rr_sci, marker='*', color='blue', s=120, label='sources_science')
    #ax.set_thetamin(-10)
    #ax.set_thetamax(100)

    # LGS
    th_HO = np.array(eval(parser['sources_HO']['Azimuth']))
    rr_HO = np.array(eval(parser['sources_HO']['Zenith']))
    th_HO = th_HO/180*np.pi
    ax.scatter(th_HO, rr_HO, marker='*', color='green', s=120, label='sources_HO')

    # NGS
    th_LO = np.array(eval(parser['sources_LO']['Azimuth']))
    rr_LO = np.array(eval(parser['sources_LO']['Zenith']))
    th_LO = th_LO/180*np.pi
    ax.scatter(th_LO, rr_LO, marker='*', color='red', s=120, label='sources_LO')

    # Set ticks position
    max_pos = np.max([rr_sci.max(), rr_HO.max(), rr_LO.max()]) + 2*ticks_interval
    ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(0,max_pos,ticks_interval)))
    r_labels = [item.get_text() for item in ax.get_yticklabels()]
    for i in range(len(r_labels)):
        if i % 2: r_labels[i]=''
    ax.set_yticklabels(r_labels, verticalalignment = "top")

    # Put weights next to science sources
    if labels is not None:
        for i,lab in enumerate(labels):
            ax.text(th_sci[i],rr_sci[i],str(lab),color='black',fontsize=11)

    # Legend
    ax.legend(loc='lower left', bbox_to_anchor=(1, 0))
    plt.show()

    return fig, ax
