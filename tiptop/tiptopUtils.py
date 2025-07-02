# some of these imports are not used here, but are used in other modules
import os
import itertools
from datetime import datetime
from collections import defaultdict
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
    
    :param hdr: FITS header object
    :type hdr: astropy.io.fits.Header
    :param key_primary: Primary key for the header keyword
    :type key_primary: str
    :param key_secondary: Secondary key for the header keyword
    :type key_secondary: str
    :param val: Value to be added to the header keyword
    :type val: str, int, float, or list
    :param iii: Optional index for the primary key
    :type iii: int, optional
    :param jjj: Optional index for the secondary key
    :type jjj: int, optional
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
    Conversion of a fits file header into a dictionary.
    Reconstructs vectors from indexed keys and handles split values.

    :param hdr: FITS header object
    :type hdr: astropy.io.fits.Header

    :return: Nested dictionary with sections and parameters
    :rtype: dict
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

            # Fix 1: Handle key extraction properly to avoid truncation
            if ext_indx == -1:
                # No '+' found, take the rest of the string
                curr_key += key[separator_index+1:]
            else:
                # '+' found, take up to the '+'
                curr_key += key[separator_index+1:ext_indx]

            val_str = str(hdr[key])
            val_last_index = val_str.find(APPEND_TOKEN)

            # Handle split values
            if val_last_index != -1:
                curr_value += val_str[:val_last_index]
            else:
                curr_value += val_str

            # Store value when no '+' (final part)
            if ext_indx == -1:
                my_data_map[section].update({curr_key.strip(): curr_value.strip()})
                curr_value = ''
                curr_key = ''

    # Fix 2: Reconstruct vectors from indexed parameters
    result = {}
    for section, params in my_data_map.items():
        result[section] = {}

        # Group parameters by base name to reconstruct vectors
        vector_groups = defaultdict(dict)

        for param_name, value in params.items():
            # Check if parameter ends with a number (index)
            parts = param_name.split()
            if len(parts) > 1 and parts[-1].isdigit():
                # It's an indexed parameter like "Cn2Heights 0"
                base_name = ' '.join(parts[:-1])
                index = int(parts[-1])
                vector_groups[base_name][index] = value
            else:
                # Regular scalar parameter
                result[section][param_name] = value

        # Convert indexed parameters to ordered lists
        for base_name, indexed_values in vector_groups.items():
            if indexed_values:
                # Find the maximum index to create the list
                max_index = max(indexed_values.keys())
                # Create ordered list with None for missing indices
                vector = []
                for i in range(max_index + 1):
                    vector.append(indexed_values.get(i, None))
                result[section][base_name] = vector

    return result

def plot_directions(parser, ticks_interval=5, labels=None, LO_labels=None,
                    science=True, max_pos=None, add_legend=True):
    '''
    Polar plot with science and GS (sources_HO and sources_LO) directions

    :param parser: required, parameters object
    :type parser: configparser object
    :param ticks_interval: optional default=5, size of interval for ticks in the figure 
    :type ticks_interval: int
    :param labels: optional default=None, list of strings to be plotted next to the science sources
    :type labels: list
    :param LO_labels: optional default=None, list of strings to be plotted next to the LO sources
    :type LO_labels: list
    :param science: optional default=True, activate plot of science sources
    :type science: bool
    :param max_pos: optional default=None, maximum distance from axis
    :type max_pos: float
    :param add_legend: optional default=True, activate legend
    :type add_legend: bool
    
    :return: fig, ax
    :rtype: objects
    '''

    fig = plt.figure('SOURCES DIRECTIONS', figsize=(6,6))
    ax = fig.add_subplot(111, polar=True)
    ax.tick_params(labelsize=10)
    ax.set_rlabel_position(225) # theta position of the radius labels
    
    # SCIENCE
    if science:
        if 'PSF_DIRECTIONS' in parser.sections(): # For retro-compatibility
            th_sci = np.array(eval(parser['PSF_DIRECTIONS']['ScienceAzimuth']))
            rr_sci = np.array(eval(parser['PSF_DIRECTIONS']['ScienceZenith']))
        else:
            th_sci = np.array(eval(parser['sources_science']['Azimuth']))
            rr_sci = np.array(eval(parser['sources_science']['Zenith']))
        th_sci = th_sci/180*np.pi
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
    if max_pos is None:
        if science:
            max_pos = np.max([rr_sci.max(), rr_HO.max(), rr_LO.max()]) + 2*ticks_interval
        else:
            max_pos = np.max([rr_HO.max(), rr_LO.max()]) + 2*ticks_interval
    ax.set_ylim(0, max_pos)
    ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(0,max_pos,ticks_interval)))
    r_labels = [item.get_text() for item in ax.get_yticklabels()]
    for i in range(len(r_labels)):
        if i % 2: r_labels[i]=''
    ax.set_yticklabels(r_labels, verticalalignment = "top")

    # Put weights next to science sources
    if labels is not None:
        for i,lab in enumerate(labels):
            ax.text(th_sci[i],rr_sci[i],str(lab),color='black',fontsize=11)
    if LO_labels is not None:
        for i,lab in enumerate(LO_labels):
            ax.text(th_LO[i],rr_LO[i],str(lab),color='black',fontsize=11)

    # Legend
    if add_legend:
        ax.legend(loc='lower left', bbox_to_anchor=(1, 0))
    plt.show()

    return fig, ax
