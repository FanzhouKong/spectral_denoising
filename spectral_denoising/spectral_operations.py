import re
import pandas as pd
import itertools
from tqdm import tqdm
import numpy as np
import scipy
import os
import ms_entropy as me
import math

_numpy_formula_format = np.int16
def normalized_entropy(msms):
    return (spectral_entropy(msms)/math.log(len(msms)))**4
def spectral_entropy(msms, pmz = None):
    """
    Calculate the entropy of the givens.

    Parameters:
        msms (numpy.ndarray): A 2D array where the each item are formated as [pmz, intensity].

        pmz (float, optional): The precursor m/z value. If provided, precursors in both spectra will be removed.
    Returns:
        float: The entropy of the query MS/MS spectrum.
    """
    if isinstance(msms, float):
        return np.nan
    if pmz is not None:
        msms = truncate_spectrum(msms, pmz-1.6)
    S = me.calculate_spectral_entropy(msms)
    return S
def entropy_similairty(msms1, msms2,pmz=None, ms2_error = 0.02):
    """
    Calculate the entropy similarity between two mass spectrometry spectra.

    Parameters:
        msms1 (numpy.ndarray): The first mass spectrometry spectrum. If a float is provided, NaN is returned.
        msms1 (numpy.ndarray): The second mass spectrometry spectrum. If a float is provided, NaN is returned.
        pmz (float, optional): The precursor m/z value. If provided, precursors in both spectra will be removed.
        ms2_error (float, optional): The tolerance for matching peaks in the spectra. Default is 0.02.
    Returns:
        float: The entropy similarity between the two spectra. Returns NaN if either input spectrum is invalid.
    """

    if isinstance(msms1, float) or isinstance(msms2, float):
        return np.nan
    if pmz is not None:
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
    """
    Compare two mass spectra and return the spectrum of the second input 
    that does not overlap with the first input. Juist a helper function, not actually in use.

    Args:
        msms1 (numpy.ndarray): The first mass spectrum to compare.
        msms2 (numpy.ndarray): The second mass spectrum to compare.
    Returns:
        numpy.ndarray: A packed spectrum of mass and intensity values from `msms2` 
              that do not overlap with `msms1`.
    """
    if len(msms2)<len(msms1):
        msms_temp = msms2
        msms2 = msms1
        msms1 = msms_temp
    mass1, intensity1 = break_spectrum(msms1)
    mass2, intensity2 = break_spectrum(msms2)

    indices = [index for index, item in enumerate(mass2) if item not in mass1]
    return pack_spectrum(mass2[indices], intensity2[indices])
def search_ions(msms, mz, span = 3):
    """
    Search for ions within a specified mass-to-charge ratio (m/z) range in a given mass spectrum.

    Parameters:
        msms (numpy.ndarray): The mass spectrum data.
        mz (float): The target mass-to-charge ratio to search for.
        span (float, optional): The range around the target m/z to search within. Default is 3.
    Returns:
        numpy.ndarray: the slice of MS/MS spectra in the given region.
    """

    mass, intensity = break_spectrum(msms)
    idx_left, idx_right = mass.searchsorted([mz-span, mz+span])
    return(pack_spectrum(mass[idx_left:idx_right], intensity[idx_left:idx_right]))

def break_spectrum(spectra):
    """
    Breaks down a given spectrum into its mass and intensity components. Not often used.

    Parameters:
        spectra (numpy.ndarray): The input spectrum data. If a np.nan is provided, it returns two empty lists.
    Returns:
        numpy.ndarray: A MS/MS spectrum formated in 2D array where the each item are formated as [pmz, intensity].
    """

    if isinstance(spectra, float):
        return ([],[])
    spectra = np.array(spectra)
    mass = spectra.T[0]
    intensity = spectra.T[1]
    return mass, intensity

def pack_spectrum(mass, intensity):
    """
    Inverse of break_spectrum. Packs mass and intensity arrays into a single 2D array, which is standardized MS/MS spectrum data format in this project.
    This function takes two arrays, `mass` and `intensity`, and combines them into a single 2D array where each row 
    corresponds to a pair of mass and intensity values. If either of the input arrays is empty, the function returns NaN.
    
    Parameters:
        mass (numpy.ndarray): An array of mass values.
        intensity (numpy.ndarray): An array of intensity values.
    Returns:
        numpy.ndarray: A 2D array with mass and intensity pairs if both input arrays are non-empty, otherwise NaN.
    """

    if len(mass)>0 and len(intensity)>0:
        return(np.array([mass, intensity]).T)
    else:
        return(np.nan)
def add_spectra(msms1, msms2):
    """
    Add two spectra together.
    This function takes two spectra (msms1 and msms2) and combines them. If one of the inputs is a float and the other 
    is not, it returns the non-float input. If both inputs are floats, it returns NaN.

    Parameters:
        msms1 (numpy.ndarray): The first spectrum.
        msms2 (numpy.ndarray): The second spectrum.
    Returns:
        numpy.ndarray: The combined spectrum if both inputs are not floats, one of the inputs if the other is a float, or NaN if both inputs are floats.
    Notes:
        - This function is very naive mixing of 2 spectrum. If you wished to formulate the intensity, please do it before using this function.
                    
    """

    if isinstance(msms1, float) == False and isinstance(msms2, float) == False:
        
        return(sort_spectrum(np.concatenate([msms1, msms2])))
    if isinstance(msms1, float) and isinstance(msms2, float) == False:
        return msms2
    elif isinstance(msms2, float) and isinstance(msms1, float) == False:
        return msms1
    else:
        return np.nan
def normalize_spectrum(msms):
    """
    Normalize the intensity values of a given mass spectrum.
    This function takes a mass spectrum (msms) as input, transposes it, and normalizes
    the intensity values (second row) by dividing each intensity by the sum of all intensities.
    The normalized spectrum is then transposed back to its original form and returned.

    Parameters:
        msms (numpy.ndarray): A 2D numpy array where the first row contains mass-to-charge ratios (m/z) and the second row contains intensity values.
    Returns:
        numpy.ndarray: A 2D numpy array with the same shape as the input, where the intensity values have been normalized.
    """

    msms_T = msms.T
    msms_T[1]=np.array([msms_T[1][i]/np.sum(msms_T[1]) for i in range(0, len(msms_T[1]) )])
    return msms_T.T 
def sort_spectrum(msms):
    """
    Sorts the spectrum data based on m/z values.

    Parameters:
        msms (numpy.ndarray): A 2D numpy array.
    Returns:
        numpy.ndarray: A 2D numpy array with the same shape as the input, but sorted by the m/z values in ascending order.
    """
    if isinstance(msms, float) or len(msms) == 0:
        return np.nan
    msms_T = msms.T
    order = np.argsort(msms_T[0])
    msms_T[0] = msms_T[0][order]
    msms_T[1] = msms_T[1][order]

    return msms_T.T
def remove_precursor(msms, pmz = None):
    """
    Removes the precursor ion from the given mass spectrometry/mass spectrometry (MS/MS) spectrum.

    Parameters:
        msms (numpy.ndarray): A 2D numpy array.
        pmz (float, optional): The precursor m/z value. If not provided, function will try to guess from the spectrum.
    Returns:
        numpy.ndarray: The truncated MS/MS spectrum with the precursor ion removed.
    """

    if isinstance(msms, float):
        return np.nan
    if pmz is None:
        pmz = max(break_spectrum(msms)[0])
    msms_t = truncate_spectrum(msms, pmz-1.6)
    return msms_t
def sanitize_spectrum(msms):
    """
    Sanitize the given mass spectrum.
    This function performs the following operations on the input mass spectrum:
    1. If the input is a nan, it returns NaN.
    2. Sorts the spectrum using the `sort_spectrum` function.
    3. Removes zero intensity ions using the `remove_zero_ions` function.

    Parameters:
        msms (numpy.ndarray): The mass spectrum to be sanitized. 
    Returns:
        numpy.ndarray: The sanitized mass spectrum. If the input is a nan, returns nan.
    """

    if isinstance(msms, float):
        return np.nan
    msms = sort_spectrum(msms)
    msms = remove_zero_ions(msms)
    return msms
def truncate_spectrum(msms, max_mz):
    """
    Truncate the given mass spectrum to only include peaks with m/z values less than or equal to max_mz.

    Parameters:
        msms (numpy.ndarray): The mass spectrum to be truncated. If it is an empty spectrum (np.nan), will also return np.nan.
        max_mz (float): The maximum m/z value to retain in the truncated spectrum.
    Returns:
        numpy.ndarray: The truncated mass spectrum with m/z values less than or equal to max_mz.

    """


    if isinstance(msms, float):
        return np.nan
    msms = sort_spectrum(msms)
    mass, intensity = msms.T[0], msms.T[1]
    upper_allowed=np.searchsorted(mass, max_mz,side = 'left')
    mass = mass[:upper_allowed]
    intensity = intensity[:upper_allowed]
    return pack_spectrum(mass, intensity)

def slice_spectrum(msms, break_mz):
    """
    Slices a mass spectrum into two parts based on a given m/z value.

    Parameters:
        msms (numpy.ndarray): The mass spectrum data, where each row represents a peak with m/z and intensity values. If a empty spectrum is provided, the function returns NaN.
        break_mz (float): The break point where to slice the spectrum.
    Returns:
        tuple: A tuple containing two numpy.ndarrays:
            - The first array contains all peaks with m/z values less than the break_mz.
            - The second array contains all peaks with m/z values greater than or equal to the break_mz.
    """

    if isinstance(msms, float):
        return np.nan
    
    idx = np.searchsorted(msms.T[0], break_mz, side = 'left')
    return(msms[:idx], msms[idx:])
def standardize_spectrum(ms):
    """
    Standardizes the intensity values of a given mass spectrum so that the base peak will have intensity of 1.

    Parameters:
        ms (numpy.ndarray): A 2D array where the first column represents mass values and the second column represents intensity values.
    Returns:
        numpy.ndarray: A 2D array with the same mass values and standardized intensity values. The intensity values are normalized to the range [0, 1] and rounded to 4 decimal places.
    """
    
    mass, intensity = ms.T[0],ms.T[1]
    intensity = intensity/np.max(intensity)
    intensity = np.round(intensity, 4)
    return(np.array([mass, intensity]).T)
    
def remove_zero_ions(msms):
    """
    Remove zero intensity ions from a mass spectrometry dataset.

    Parameters:
        msms (numpy.ndarray or float): MS/MS spectrum in 2D numpy array.
    Returns:
        numpy.ndarray: A filtered 2D numpy array with rows where the second column (ion intensities) is greater than zero, or np.nan if the input is an empty spectrum.
    """

    if isinstance(msms, float) or len(msms) == 0:
        return np.nan
    to_keep = msms.T[1] > 0
    return msms[to_keep]

def arr_to_str(msms):
    '''
    helper function for read_df and save_df
    '''
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
    '''
    helper function for read_df and save_df
    '''
    if isinstance(msms, float):
        return np.nan
    spec_raw = np.array([x.split('\t') for x in msms.split('\n')], dtype=np.float32)
    return(spec_raw)









