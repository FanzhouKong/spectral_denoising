#!/usr/bin/env python
# coding: utf-8

# In[ ]:
__all__ = []

import re
import pandas as pd
import itertools
from tqdm import tqdm
import numpy as np
import scipy
import os
import ms_entropy as me
import bisect
import warnings
import math
warnings.filterwarnings("ignore")
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import warnings
warnings.filterwarnings("ignore")
# import toolsets.denoising_related_functions as de
# from toolsets.search import quick_search_values
_numpy_formula_format = np.int16
def spctrum_entropy(msms):
    intensity = msms.T[1]
    S = scipy.stats.entropy(intensity)
    return S
def entropy_similairty(msms1, msms2,pmz=None, ms2_error = 0.02):
    if isinstance(msms1, float) or isinstance(msms2, float):
        return np.nan
    msms1 = truncate_spectrum(msms1, pmz-1.6)
    msms2 = truncate_spectrum(msms2, pmz-1.6)
    if isinstance(msms1, float) or isinstance(msms2, float):
        return np.nan
    if pmz is not None:
        similarity = me.calculate_entropy_similarity(msms1, msms2, ms2_tolerance_in_da = ms2_error, noise_threshold=0.00, clean_spectra=True,max_mz = pmz-1.6)
    else:
        similarity = me.calculate_entropy_similarity(msms1, msms2, ms2_tolerance_in_da = ms2_error, noise_threshold=0.00, clean_spectra=True)
    return similarity
def compare_spectra(msms1, msms2):
    if len(msms2)<len(msms1):
        msms_temp = msms2
        msms2 = msms1
        msms1 = msms_temp
    mass1, intensity1 = break_spectrum(msms1)
    mass2, intensity2 = break_spectrum(msms2)

    indices = [index for index, item in enumerate(mass2) if item not in mass1]
    return pack_spectrum(mass2[indices], intensity2[indices])
def search_ions(msms, mz, span = 3):
    mass, intensity = break_spectrum(msms)
    idx_left, idx_right = mass.searchsorted([mz-span, mz+span])
    return(pack_spectrum(mass[idx_left:idx_right], intensity[idx_left:idx_right]))

def break_spectrum(spectra):
    if isinstance(spectra, float):
        return ([],[])
    spectra = np.array(spectra)
    mass = spectra.T[0]
    intensity = spectra.T[1]
    return mass, intensity

def pack_spectrum(mass, intensity):
    if len(mass)>0 and len(intensity)>0:
        return(np.array([mass, intensity]).T)
    else:
        return(np.nan)
def add_spectra(msms1, msms2):
    if isinstance(msms1, float) == False and isinstance(msms2, float) == False:
        
        return((np.concatenate([msms1, msms2])))
    if isinstance(msms1, float) and isinstance(msms2, float) == False:
        return msms2
    elif isinstance(msms2, float) and isinstance(msms1, float) == False:
        return msms1
    else:
        return np.nan
def normalize_spectrum(msms):
    msms_T = msms.T
    msms_T[1]=np.array([msms_T[1][i]/np.sum(msms_T[1]) for i in range(0, len(msms_T[1]) )])
    return msms_T.T 
def sort_spectrum(msms):
    msms_T = msms.T
    order = np.argsort(msms_T[0])
    msms_T[0] = msms_T[0][order]
    msms_T[1] = msms_T[1][order]

    return msms_T.T
def remove_precursor(msms, pmz = None):
    if isinstance(msms, float):
        return np.nan
    if pmz is None:
        pmz = max(break_spectrum(msms)[0])
    msms_t = truncate_spectrum(msms, pmz-1.6)
    return msms_t
def sanitize_spectrum(msms):
    if isinstance(msms, float):
        return np.nan
    msms = sort_spectrum(msms)
    msms = remove_zero_ions(msms)
    return msms
def truncate_spectrum(msms, max_mz):
    if isinstance(msms, float):
        return np.nan
    msms = sort_spectrum(msms)
    mass, intensity = msms.T[0], msms.T[1]
    upper_allowed=np.searchsorted(mass, max_mz,side = 'left')
    mass = mass[:upper_allowed]
    intensity = intensity[:upper_allowed]
    return pack_spectrum(mass, intensity)

def slice_spectrum(msms, break_mz):
    if isinstance(msms, float):
        return np.nan
    
    idx = np.searchsorted(msms.T[0], break_mz, side = 'left')
    return(msms[:idx], msms[idx:])
def standardize_spectrum(ms):
    mass, intensity = ms.T[0],ms.T[1]
    intensity = intensity/np.max(intensity)
    intensity = np.round(intensity, 4)
    return(np.array([mass, intensity]).T)
    
def remove_zero_ions(msms):
    if isinstance(msms, float):
        return np.nan
    to_keep = msms.T[1] > 0
    return msms[to_keep]






def arr_to_str(msms):
    if isinstance(msms, float):
        return np.nan
    mass = []
    intensity = []
    for n in range(0, len(msms)):
        mass.append(msms[n][0])
        intensity.append(msms[n][1])
    if len(mass)>0 and len(intensity)>0 and len(mass)==len(intensity):
        intensity_return = [str(inten) + '\n' for (inten) in (intensity[:-1])]
        intensity_return.append(str(intensity[-1]))
        mass_cali_tab = [str(mas) + '\t' for (mas) in mass]
        list_temp = [None]*(len(mass_cali_tab)+len(intensity_return))
        list_temp[::2] = mass_cali_tab
        list_temp[1::2] = intensity_return
        list_temp = ''.join(list_temp)
        return(list_temp)
    else:
        return(np.nan)
def str_to_arr(msms):
    if isinstance(msms, float):
        return np.nan
    spec_raw = np.array([x.split('\t') for x in msms.split('\n')], dtype=np.float32)
    return(spec_raw)









