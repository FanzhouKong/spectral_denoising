#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import re
import pandas as pd
import spectral_entropy as se
import itertools
from tqdm import tqdm
import numpy as np
import scipy.stats
import os
import bisect
import warnings
import math
warnings.filterwarnings("ignore")
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import toolsets.denoising_related_functions as de
from toolsets.search import quick_search_values

def remove_zero_ions(ms1):
    mass1, intensity1 = break_spectra(ms1)
    intensity_array = np.array(intensity1)
    mass_array = np.array(mass1)
    to_keep = intensity_array > 0
    mass_1_cleaned = list(mass_array[to_keep])
    intensity_1_cleaned = list(intensity_array[to_keep])
    ms1_cleaned = pack_spectra(mass_1_cleaned, intensity_1_cleaned)
    return ms1_cleaned
def check_spectrum(msms):
    try:
        mass, intensity = break_spectra(msms)
    except:
        print("it is not a spectrum!")
        return()
    if len(mass)<1 or len(intensity)<1 or max(intensity)==0:
        print("spectrum is empty!")
        return(-1)
    return(sum(intensity))

def sort_spectrum(msms):
    if type(msms) == str:
        msms = ast.literal_eval(msms)

    mass, intensity = break_spectra(msms)
    order = np.argsort(mass)
    mass = mass[order]
    intensity = intensity[order]
    return pack_spectra(mass, intensity)


# def clean_spectrum(msms,
#                     max_mz: float = None,
#                     keep_precursor = False,
#                     tolerance:float = 0.02,
#                     ifppm: bool = False,
#                     noise_level: float = 0.000):
#     # remove precursor peak
#     if msms !=msms:
#         return (np.NaN)
#     msms = sort_spectrum(msms)
#     max_mz = float(max_mz)
#     tol = set_tolerance(mass_error=tolerance,
#         ifppm = ifppm,precursormz = max_mz)
#     if tolerance is None:
#         raise RuntimeError("MS2 tolerance need to be set!")
#     if max_mz is not None:
#         if keep_precursor == False:
#             msms = truncate_msms(msms, max_mz-1.5)
#         else:
#             msms = truncate_msms(msms, max_mz)
#     # bin spectra
#     if msms is not np.NAN:
#         msms = bin_spectrum(msms=msms, precursormz = max_mz, tol = tol)
#         # remove noise by level
#         # msms = denoising_by_threshold(msms=msms, threshold = noise_level)
#         msms = sort_spectrum(msms)
#         msms = normalize_spectrum(msms)
#         return(msms)
#     else:
#         return(msms)
def get_fragment_intensity(msms, pmz):

    msms_frag = truncate_spectrum(msms, pmz-1.6)
    if isinstance(msms_frag, float):
        return 0
    mass, intensity = break_spectra(msms_frag)
    return np.sum(intensity)
def truncate_spectrum(msms, max_mz):
    if isinstance(msms, float):
        return np.NAN
    msms = sort_spectrum(msms)
    mass, intensity = break_spectra(msms)
    upper_allowed=np.searchsorted(mass, max_mz,side = 'left')
    # upper_allowed = bisect.bisect_left(mass, max_mz)
    # pep_index = mass.index(parention)
    mass_frag = mass[0:upper_allowed]
    intensity_frag = intensity[0:upper_allowed]
    if len(mass_frag)>0 and len(intensity_frag)>0:
        return((pack_spectra(mass_frag, intensity_frag)))
    else:
        return(np.NAN)


def add_spectra(msms1, msms2):
    mass1,intensity1 = break_spectra(msms1)
    mass2, intensity2 = break_spectra(msms2)
    mass1.extend(mass2)
    intensity1.extend(intensity2)
    msms_mix = pack_spectra(mass1, intensity1)
    msms_mix = sort_spectrum(msms_mix)
    # msms_mix_binned = bin_spectrum(msms_mix)
    return(msms_mix)
def bin_spectrum(msms,mass_error = 0.01):
    # normalize spectrum is included already!!!
    mass, intensity = break_spectra(sort_spectrum(msms))
    mass_binned = []
    intensity_binned = []
    while len(mass)>0:
        idx = np.argmax(intensity)
        idx_start = np.searchsorted(mass, mass[idx]-mass_error, 'left')
        idx_end = np.searchsorted(mass, mass[idx]+mass_error, 'right')
        mass_temp = mass[idx_start:idx_end]
        intensity_temp = intensity[idx_start:idx_end]
        mass_binned.append(mass[idx])
        intensity_binned.append(np.sum(intensity_temp))
        mass = np.concatenate((mass[:idx_start], mass[idx_end:]))
        intensity = np.concatenate((intensity[:idx_start], intensity[idx_end:]))
    msms_binned = sort_spectrum(pack_spectra(mass_binned, intensity_binned))
    return msms_binned
def clean_spectrum(msms):
    msms = bin_spectrum(msms)
    msms = remove_zero_ions(msms)
    msms = normalize_spectrum(msms)
    msms = sort_spectrum(msms)
    return(msms)

# mass_temp, intensity_temp = break_spectra(msms)
    # bin_left = pd.DataFrame({'mass': mass_temp, 'intensity': intensity_temp})
    # mass_bin = []
    # intensity_bin = []
    # while(len(bin_left)>0):
    #
    #     binn, bin_left = make_bin(bin_left, tol)
    #
    #     # mass_bin.append(round(bin['mass'].loc[bin['intensity']== bin['intensity'].max()], 6))
    #     # max_index = bin['intensity'].idxmax()
    #     mass_bin.append(round(binn['mass'].median(), 6))
    #     # mass_bin.append(round(bin['mass'][max_index], 6))
    #     # intensity_bin.append(round(bin['intensity'].max(),6))
    #     intensity_bin.append(round(binn['intensity'].sum(),6))
    # msms_bin = normalize_spectrum(sort_spectrum(pack_spectra(mass_bin, intensity_bin)))
    #
    # return(msms_bin)
def make_bin(bin_left, tol, if_sorted = True):
    if if_sorted == False:
        bin_left.sort_values(by='mass', ascending=True, inplace=True)
    # bin_left needed to be sorted!!!
    max_index = bin_left['intensity'].idxmax()
    # step = tol/2
    step = tol
    binn =quick_search_values(bin_left, 'mass', bin_left['mass'][max_index]-step, bin_left['mass'][max_index]+step)
    bin_left_return = bin_left.drop(binn.index)
    return(binn, bin_left_return)
def weighted_average_spectra(data_subset, typeofmsms = 'peaks', weight_col = 'ms1_precursor_intensity'):
    # if(len(data_subset)<2):
    #     print("you cannot make weighted average spectra using only 1 spectra")
    #     return(np.NAN)
    # precursormz = float(data_subset.iloc[0]['reference_precursor_mz'])
    # msms_com = []

    mass_com = []
    intensity_com = []
    ms1_intensity =[]
    sum = 0
    for index, row in data_subset.iterrows():
        if (row[typeofmsms])==(row[typeofmsms]) and row[weight_col]!=0:
            mass_temp, intensity_temp = break_spectra(normalize_spectrum(row[typeofmsms]))
        else:
            continue
        sum = sum+row[weight_col]
        # msms_bin_temp =normalize_spectra(row[typeofmsms])
        # msms_bin_temp = bin_spectrum(row[typeofmsms])

        mass_com.extend(mass_temp)
        intensity_com.extend(intensity_temp)
        # num_peaks_com.extend(len(mass_temp))
        ms1_intensity.extend([row[weight_col]]*len(mass_temp))# five places to change
        # sum = sum +row['intensity']
    if len(mass_com)==0:
        return(np.NAN)

    bin_left = pd.DataFrame({'mass': mass_com, 'intensity': intensity_com, 'ms1_precursor_intensity':ms1_intensity})
    bin_left['intensity_weighted']=[x*y/sum for x, y in zip(bin_left['intensity'], bin_left['ms1_precursor_intensity'])]
    # return(bin_left)
    bin_left.sort_values(by = 'mass', ascending=True, inplace=True)
    msms_weighted = bin_spectrum(pack_spectra(bin_left['mass'].tolist(), bin_left['intensity_weighted'].tolist()))
    msms_weighted = normalize_spectrum(msms_weighted)
    return msms_weighted

# def weighted_average_spectra(data_subset, typeofmsms = 'peaks', mass_error = 0.02, weight_col = 'ms1_precursor_intensity'):
#     # if(len(data_subset)<2):
#     #     print("you cannot make weighted average spectra using only 1 spectra")
#     #     return(np.NAN)
#     # precursormz = float(data_subset.iloc[0]['reference_precursor_mz'])
#     # msms_com = []
#     mass_com = []
#     intensity_com = []
#     ms1_intensity =[]
#     sum = 0
#     for index, row in data_subset.iterrows():
#         if (row[typeofmsms])==(row[typeofmsms]) and row[weight_col]!=0:
#             mass_temp, intensity_temp = break_spectra(row[typeofmsms])
#         else:
#             continue
#         sum = sum+row[weight_col]
#         # msms_bin_temp =normalize_spectra(row[typeofmsms])
#         # msms_bin_temp = bin_spectrum(row[typeofmsms])
#
#         mass_com.extend(mass_temp)
#         intensity_com.extend(intensity_temp)
#         # num_peaks_com.extend(len(mass_temp))
#         ms1_intensity.extend([row[weight_col]]*len(mass_temp))# five places to change
#         # sum = sum +row['intensity']
#     if len(mass_com)==0:
#         return(np.NAN)
#     bin_left = pd.DataFrame({'mass': mass_com, 'intensity': intensity_com, 'ms1_precursor_intensity':ms1_intensity})
#     bin_left.sort_values(by='mass', ascending=True, inplace=True)
#     # return(bin_left)
#     # bin_left['adjusted_intensity']=bin_left['intensity']*bin_left['ms1_precursor_intensity_ratio']
#
#     # return(bin_left)
#     # msms_binned = pack_spectra(bin_left['mass'].tolist(), bin_left['adjusted_intensity'].tolist())
#     # msms_binned = normalize_spectra(msms_binned)
#     # return(bin_left)
#     mass_consensus = []
#     intensity_consensus =[]
#     # return(bin_left)
#     # print('i am in new method')
#     while(len(bin_left)>0):
#         binn, bin_left = make_bin(bin_left, mass_error)
#         mass_temp = 0
#         intensity_temp = 0
#         denominator = binn['ms1_precursor_intensity'].sum()
#         for index, row in binn.iterrows():
#             mass_temp = mass_temp+row['mass']*row['ms1_precursor_intensity']
#             intensity_temp = intensity_temp+row['intensity']*row['ms1_precursor_intensity']
#         mass_consensus.append(mass_temp/denominator)
#         intensity_consensus.append(intensity_temp/denominator)
#
#
#         # if len(binn)>1:
#         #     temp_mass = (binn['mass']*binn['ms1_precursor_intensity_ratio']).sum()
#         #     temp_intensity =(binn['intensity']*binn['ms1_precursor_intensity_ratio']).sum()
#         # else:
#         #     temp_mass = binn['mass'].sum()
#         #     temp_intensity = (binn['intensity']*binn['ms1_precursor_intensity_ratio']).sum()
#         # mass_consensus.append(round(temp_mass,6))
#         # intensity_consensus.append(round(temp_intensity,6))
#
#     msms_consensus = sort_spectrum(pack_spectra(mass_consensus, intensity_consensus))
#
#     msms_consensus = bin_spectrum(msms_consensus, tol = 0.02)
#     return(msms_consensus)










def denoising_by_threshold(msms, threshold = 0.005):
    mass_raw, intensity_raw = break_spectra(msms)
    threshold_use = max(intensity_raw)*threshold
    idx=[index for (index, number) in enumerate(intensity_raw) if number > threshold_use]
    intensity_updated = [intensity_raw[i] for i in idx]
    mass_updated = [mass_raw[i] for i in idx]
    msms_updated = pack_spectra(mass_updated, intensity_updated)
    return(msms_updated)

def guess_charge(adduct):
    if adduct[-1]=="+":
        return('pos')
    elif adduct[-1]=='-1':
        return('neg')
    else:
        return(np.NAN)

def standardize_spectra(ms):
    mass, intensity = break_spectra(ms)
    intensity = intensity/np.max(intensity)
    intensity = np.round(intensity, 4)
    return(pack_spectra(mass, intensity))
    
def break_spectra(spectra):
    if isinstance(spectra, float):
        return ([],[])
    spectra = np.array(spectra)
    mass = spectra.T[0]
    intensity = spectra.T[1]
    return mass, intensity

def pack_spectra(mass, intensity):
    if len(mass)>0 and len(intensity)>0:
        return(np.array([[mass[i], intensity[i]] for i in range(0, len(mass))], dtype=np.float32))
    else:
        return(np.NAN)


def convert_arr_to_string(msms):
    if isinstance(msms, float):
        return np.NAN
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
        return(np.NAN)
def standardize_msdial(msms):
    spec_raw = np.array([x.split(':') for x in msms.split(' ')], dtype=np.float32)
    return spec_raw
    # return(msdial_msms.replace(' ', '\n').replace(':', '\t'))
def convert_scc_to_string(msms):
    mass = []
    intensity = []
    lst = msms.split(';')
    for l in range(len(lst)):
        list_temp = lst[l].split(':')
        mass.append(float(list_temp[0]))
        intensity.append(float(list_temp[1]))
    msms_return = pack_spectra(mass, intensity)
    return(msms_return)

def convert_string_to_arr(msms):
    if isinstance(msms, float):
        return np.NAN
    spec_raw = np.array([x.split('\t') for x in msms.split('\n')], dtype=np.float32)
    return(spec_raw)

def spectral_entropy(msms):
    if isinstance(msms, float):
        return(-1)
    mass, intensity = break_spectra(msms)
    return(scipy.stats.entropy(intensity))
def normalized_entropy(msms, order = 4):
    if isinstance(msms, float):
        return(pack_spectra([],[]))
    npeak = num_peaks(msms)
    mass, intensity = break_spectra(msms)
    normalized_entropy = (scipy.stats.entropy(intensity)/math.log(npeak))**order
    if normalized_entropy==normalized_entropy:
        return (normalized_entropy)
    else:
        return(np.NAN)

def convert_lcb_to_arr(msms):
    msms_list = msms.split(' ')
    mass = []
    intensity = []
    for mi in msms_list:
        mass.append(float(mi.split(':')[0]))
        intensity.append(float(mi.split(':')[1]))
    msms_string = pack_spectra(mass, intensity)
    return(msms_string)
# def add_spectra(msms1, msms2):
#     mass1, intensity1 = break_spectra(msms1)
#     mass2, intensity2 = break_spectra(msms2)
#     mass_mix =
# def make_composite_spectra(data_subset, typeofmsms = 'msms', tolerance = 0.01, ifnormalize = False, ifppm = False):
#     precursormz = float(data_subset.iloc[0]['PRECURSORMZ'])
#     if ifppm:
#         if precursormz <400:
#             tol = 0.004
#         else:
#             tol = precursormz*(tolerance/1E6)
#     else:
#         tol = tolerance
#     mass_com = []
#     intensity_com = []
#     for index, row in data_subset.iterrows():

#         mass_temp, intensity_temp = break_spectra(normalize_spectra(row[typeofmsms]) )
#         mass_com.extend(mass_temp)
#         intensity_com.extend(intensity_temp)
#     msms_com = sort_spectra(pack_spectra(mass_com, intensity_com))
#     msms_com = bin_spectra(msms_com,precursormz, tol, ifppm =ifppm, ifnormalize=ifnormalize)
#     return((msms_com))




# def straight_adding_spectra(data_subset, typeofmsms = 'peaks'):
#     mass_com = []
#     intensity_com = []
#     for index, row in data_subset.iterrows():
#         # msms_bin_temp =normalize_spectra(row[typeofmsms])
#         # msms_bin_temp = bin_spectra(row[typeofmsms],precursormz, tol,ifnormalize=False, ifppm =ifppm)
#         mass_temp, intensity_temp = break_spectra(row[typeofmsms])
#         mass_com.extend(mass_temp)
#         intensity_com.extend(intensity_temp)
#         # num_peaks_com.extend(len(mass_temp))
#     msms_added = sort_spectra(pack_spectra(mass_com, intensity_com))
#     msms_added = normalize_spectra(msms_added)
#     return(msms_added)






# def make_consensus_spectra(data_subset, typeofmsms = 'msms', tolerance = 0.01, ifnormalize = False, ifppm = False):
#     # if(len(data_subset)<2):
#     #     print("you cannot make weighted average spectra using only 1 spectra")
#     #     return(np.NAN)
#     precursormz = float(data_subset.iloc[0]['PRECURSORMZ'])
#
#     mass_consensus = []
#     intensity_consensus =[]
#     mass_com =[]
#     intensity_com=[]
#     if ifppm:
#         if float(precursormz) <400:
#             tol = 0.004
#         else:
#             tol = float(precursormz)*(tolerance/1E6)
#     else:
#         tol = tolerance
#     # for msms in msms_s:
#     #     msms_bin_temp = bin_spectra(msms,precursormz, tolerance=tol,ifnormalize=ifnormalize, ppm =ppm)
#     #     mass_temp, intensity_temp = break_spectra(msms_bin_temp)
#     #     mass_com.extend(mass_temp)
#     #     intensity_com.extend(intensity_temp)
#     for index, row in data_subset.iterrows():
#         msms_bin_temp = bin_spectra(row[typeofmsms],precursormz, tol,ifnormalize=ifnormalize, ifppm =ifppm)
#         mass_temp, intensity_temp = break_spectra(msms_bin_temp)
#         mass_com.extend(mass_temp)
#         intensity_com.extend(intensity_temp)
#         # ms1_intensity_com.extend([row['intensity']]*len(mass_temp))
#     bin_left = pd.DataFrame({'mass': mass_com, 'intensity': intensity_com})
#
#     while(len(bin_left)>0):
#         bin, bin_left = make_bin(bin_left, tol)
#         mass_consensus.append(round(bin['mass'].median(), 6))
#         intensity_consensus.append(round(bin['intensity'].median(),6))
#     msms_consensus = sort_spectra(pack_spectra(mass_consensus, intensity_consensus))
#     msms_consensus = normalize_spectra(msms_consensus)
#     return(msms_consensus)




    # return(((scipy.stats.entropy(convert_string_to_nist(msms)[:, 1]))/math.log(npeak))^order)
# needs further modification

import ms_entropy as me
import ast
def entropy_identity(msms1, msms2,pmz, ms2_error = 0.02):
    if isinstance(msms1, float) or isinstance(msms2, float):
        return np.NAN
    if isinstance(msms1, str):
        msms1 = ast.literal_eval(msms1)
    if isinstance(msms2, str):
        msms2 = ast.literal_eval(msms2)
    similarity = me.calculate_entropy_similarity(msms1, msms2, ms2_tolerance_in_da = ms2_error, noise_threshold=0.00, clean_spectra=True,max_mz = pmz-1.6)
    return similarity

    # return(se.similarity(convert_string_to_nist(msms1), convert_string_to_nist(msms2), 'entropy',
    #     ms2_da = mass_error, need_clean_spectra = False, need_normalize_result = True))

import random
def generate_noise(pmz, lamda, n = 100):
    if int(n)!= n:
        n = np.int64(np.ceil(n))
    else:
        n = n
    # Generate a random variable from a uniform distribution in the range [a, b]
    mass = [random.uniform(50, pmz) for _ in range(n)]

    # size specifies the number of random variates to generate.

    # Generating Poisson-distributed random variables
    intensity = np.random.poisson(lam=lamda, size=n)
    intensity = intensity/100
    return(pack_spectra(mass, intensity))
def generate_chemical_noise(pmz, lamda, polarity,formula_db,n = 100):
    mass_e =  -0.00054858026
    if polarity =='+':
        coe = 1
    elif polarity =='-':
        coe = -1
    else:
        print('cannot determine adduct polarity!')
        return()
    if int(n)!= n:
        n = np.int64(np.ceil(n))
    else:
        n = n
    all_possible_mass = np.array(formula_db['mass'])
    idx_left, idx_right = all_possible_mass.searchsorted([50,pmz ])
    all_allowed_mass = all_possible_mass[idx_left:idx_right]
    if idx_right-idx_left <n:
        n = idx_right-idx_left
    # Generate a random variable from a uniform distribution in the range [a, b]
    mass = np.random.choice(all_allowed_mass, size=n, replace=False)
    mass = mass+coe*mass_e
    # mass = [random.uniform(50, pmz) for _ in range(n)]

    # size specifies the number of random variates to generate.

    # Generating Poisson-distributed random variables
    intensity = np.random.poisson(lam=lamda, size=n)
    intensity = intensity/100
    return(pack_spectra(mass, intensity))
def add_noise(msms1, noise):
    msms1 = standardize_spectra(msms1)
    mass1, intensity1= break_spectra(msms1)
    mass2, intensity2 = break_spectra(noise)
    mass = np.concatenate((mass1, mass2))
    intensity = np.concatenate((intensity1, intensity2))
    msms = pack_spectra(mass, intensity)
    # msms = bin_spectrum(msms)
    return msms





# def duplicate_handling(data, typeofmsms='msms_recalibrated',mass_error = 0.01, ifppm = False, ifnormalize = True, method = 'weighedaverage'):
#     data_unique = data.drop_duplicates(subset=['key'])
#     # print("i am in new method")
#     consensus_msms = []
#     for key in tqdm(data_unique['key'].unique()):
#         data_temp =  data.loc[data['key']==key]
#         if len(data_temp) >1:
#             if method =='weighedaverage':
#                 consensus_msms.append(weighted_average_spectra(data_temp, mass_error=mass_error, ifppm = ifppm, typeofmsms =typeofmsms))
#             # else:
#             #     consensus_msms.append(make_consensus_spectra(data_temp, tolerance=mass_error, ifppm = ifppm, typeofmsms =typeofmsms, ifnormalize=ifnormalize))
#         else:
#             consensus_msms.append(normalize_spectra(data_temp.iloc[0][typeofmsms]))


#     data_unique['msms']=consensus_msms
#     # print("i am done processing the concensus spectra")
#     return(data_unique)
import toolsets.denoising_related_functions as de

def denoising(data, typeofmsms, mass_error = 0.02, ifppm = False):
    msms_consensus_denoised = []
    comments =[]
    # print("i am in new denoising method")
    for index, row in tqdm(data.iterrows(), total = data.shape[0]):
    # break
        try:
            msms_consensus_denoised.append(de.denoise_blacklist(row, typeofmsms=typeofmsms, mass_error=mass_error))
            comments.append('denoised')
        except:
            msms_consensus_denoised.append(row[typeofmsms])
            comments.append('not denoised due to some errors')
    denoised_column = typeofmsms+"_denoised"
    data[denoised_column]=msms_consensus_denoised
    data['denoised_comments']=comments
    return(data)


# def denoising_evaluation(data, msms1 = 'peaks_recalibrated', msms2 = 'peaks_recalibrated_denoised', min_explained_intensity = 80):
#     explained_intensity = []
#     max_unassigned_intensity = []
#     for index, row in (data.iterrows()):
#         explained_intensity.append(calculate_explained_intensity(row[msms1], row[msms2]))
#         max_unassigned_intensity.append(identify_max_unassigned_intensity(row[msms1], row[msms2]))
#     data['explained_intensity']=explained_intensity
#     data['max_unassigned_intensity']=max_unassigned_intensity
#     evaluations = []
#     # print(min_explained_intensity/10)
#     for index, row in (data.iterrows()):
#         if row['explained_intensity']<min_explained_intensity/100 and row['max_unassigned_intensity']>allowed_max_unassigned_intensity:
#             evaluations.append('flagged: poor quality')
#         elif row['explained_intensity']<min_explained_intensity/100:
#             evaluations.append('flagged:low assigned intensity')
#         elif row['max_unassigned_intensity']>allowed_max_unassigned_intensity:
#             evaluations.append('flagged: high unassigned intensity')
#         else:
#             evaluations.append('good quality')
#     data['evaluations']=evaluations
#     return(data)




def num_peaks(msms):
    if msms is np.NAN:
        return(0)
    else:
        mass, intensity = break_spectra(msms)
        return(len(mass))

from operator import itemgetter

# still building
# def evaluate_spectra(msm1, msms2):
#     mass_raw, intensity_raw = break_spectra(msms1)
#     mass_dr, intensity_dr = break_spectra(msms2)
#     mass_raw_fl = [float(x) for x in mass_raw]
#     mass_dr_fl = [float(x) for x in mass_dr]
#     diff_index = [i for i, item in enumerate(mass_raw_fl) if item not in set(mass_dr_fl)]
#     mass_diff = list(itemgetter(*diff_index)(mass_raw))
#     intensity_diff = list(itemgetter(*diff_index)(intensity_raw))
#     rel_intensity_diff = [number / max([float(x) for x in intensity_raw])*100 for number in [float(y) for y in intensity_diff]]
#     rel_intensity_kept = [number / max([float(x) for x in intensity_raw])*100 for number in [float(y) for y in intensity_dr]]
#     return(mass_diff,rel_intensity_diff,rel_intensity_kept)

def normalize_spectrum(msms, total = 'unit'):
    if msms is np.NAN:
        return(msms)
    mass, intensity = break_spectra(msms)
    # mass_fl = [float(x) for x in mass]
    # intensity_fl = [float(x) for x in intensity]
    if max([float(x) for x in intensity]) == 0:
        return(np.NaN)
    if total == 'unit':
        intensity_rel = [x/np.sum(intensity) for x in intensity]
    elif total == 'half':
        intensity_rel = [x/np.sum(intensity) for x in intensity]
        intensity_rel = [x/2 for x in intensity_rel]
    # intensity_rel = [number / sum([float(x) for x in intensity]) for number in [float(y) for y in intensity]]
    # intensity_rel = [round(number, 8) for number in intensity_rel]
    return(pack_spectra(mass, intensity_rel))
def cut_msms(msms, mz_lower, mz_upper):
    if isinstance(msms, float):
        return np.NAN
    mass_temp, intensity_temp = break_spectra(msms)
    search_array=np.array(mass_temp)

    if mass_temp.all() == mass_temp.all():
        index_start, index_end = search_array.searchsorted([mz_lower, mz_upper])
        mass_return = mass_temp[index_start:index_end]
        intensity_return = intensity_temp[index_start:index_end]
        return(pack_spectra(mass_return, intensity_return))
    else:
        return(np.NAN)
def _extract_ms1_intensity(peaks, mz_lower, mz_upper):
    mass_temp, intensity_temp = break_spectra(peaks)
    if mass_temp == mass_temp:
        index_start = np.searchsorted(mass_temp, mz_lower,side = 'left')
        index_end = np.searchsorted(mass_temp, mz_upper,side = 'right')
        return(np.sum(intensity_temp[index_start: index_end]))
    else:
        return(0)
def _extract_ms1_mass(peaks, mz_lower, mz_upper):

    mass_temp, intensity_temp = break_spectra(peaks)
    if mass_temp == mass_temp:
        index_start = np.searchsorted(mass_temp, mz_lower,side = 'left')
        index_end = np.searchsorted(mass_temp, mz_upper,side = 'right')
        return(mass_temp[index_start:index_end])
    else:
        return([])
# def comparing_spectrum(msms1, msms2):
#     if(num_peaks(msms1)<num_peaks(msms2)):
#         temp_msms = msms1
#         msms1 = msms2
#         msms2 = temp_msms
#     mass_raw, intensity_raw = break_spectra(msms1)
#     mass_dr, intensity_dr = break_spectra(msms2)
#     mass_raw_fl = [float(x) for x in mass_raw]
#     mass_dr_fl = [float(x) for x in mass_dr]
#     diff_index = [i for i, item in enumerate(mass_raw_fl) if item not in set(mass_dr_fl)]
#     mass_diff = list(itemgetter(*diff_index)(mass_raw))
#     intensity_diff = list(itemgetter(*diff_index)(intensity_raw))
#     rel_intensity_diff = [number / max([float(x) for x in intensity_raw])*100 for number in [float(y) for y in intensity_diff]]
#     rel_intensity_kept = [number / max([float(x) for x in intensity_raw])*100 for number in [float(y) for y in intensity_dr]]
#     return(mass_diff,rel_intensity_diff,rel_intensity_kept)


def get_top_n(peaks, n, pmz):
    peaks = truncate_spectrum(peaks, pmz -1.6)
    if num_peaks(peaks)<n:
        return(np.NAN, np.NAN)
    mass, intensity = break_spectra(peaks)
    intensity_sorted, mass_sorted = zip(*sorted(zip(intensity, mass)))
    proportion = intensity_sorted[-n]/np.sum(intensity)
    return(mass_sorted[-n], proportion)
def search_mrm_ions(msms, mrm_ion):
    mass, intensity = break_spectra(msms)
    offset = [abs(x-mrm_ion) for x in mass]
    idx = np.argmin(offset)
    if offset[idx]>0.005:
        return(mrm_ion, 0)
    else:
        return mass[idx],intensity[idx]/np.sum(intensity)
# from toolsets.denoising_related_functions import remove_precursor
# below is some evaluation functions
def calculate_explained_intensity(msms1, msms2):
    if msms1 is np.NAN or msms2 is np.NAN:
        # print('somethign goies wrong')
        return(0)
    if(num_peaks(msms1)<num_peaks(msms2)):
        temp_msms = msms1
        msms1 = msms2
        msms2 = temp_msms

    mass_frag_raw, intensity_frag_raw = break_spectra(msms1)
    # mass_frag_raw, intensity_frag_raw= break_spectra(truncate_msms(msms1, parent_ion)) 
    mass_frag_dr, intensity_frag_dr = break_spectra(msms2)
    # mass_frag_dr, intensity_frag_dr = break_spectra(truncate_msms(msms2, parent_ion))
    if sum(intensity_frag_raw)==0:
        return(np.NaN)
    else:
        return(round(sum(intensity_frag_dr)/sum(intensity_frag_raw)*100, 5))

def identify_max_unassigned_intensity(msms1, msms2,parent_ion):
    if(num_peaks(msms1)<num_peaks(msms2)):
        temp_msms = msms1
        msms1 = msms2
        msms2 = temp_msms
    mass_raw, intensity_raw = break_spectra(msms1)
    mass_frag_raw, intensity_frag_raw, mass_precursor_raw, intensity_precursor_raw = remove_precursor(mass_raw,intensity_raw, parent_ion)
    mass_dr, intensity_dr = break_spectra(msms2)
    mass_frag_dr, intensity_frag_dr, mass_precursor_dr, intensity_precursor_dr = remove_precursor(mass_dr,intensity_dr, parent_ion)

    diff_index = [i for i, item in enumerate(mass_frag_raw) if item not in set(mass_frag_dr)]
    # return(diff_index)
    if(len(diff_index)>1):
        # print("i am in wrong if")
        intensity_diff = list(itemgetter(*diff_index)(intensity_frag_raw))
        return(max(intensity_diff))
    elif(len(diff_index)==1):
        # print("i am in right if")
        return(intensity_frag_raw[diff_index[0]])
    else:
        return(0)
    #
    # return(mass_diff)

def identify_unassigned_intensity(msms1, msms2):
    if(num_peaks(msms1)<num_peaks(msms2)):
        temp_msms = msms1
        msms1 = msms2
        msms2 = temp_msms
    mass_raw, intensity_raw = break_spectra(msms1)
    mass_de, intensity_de = break_spectra(msms2)
    diff_index = [i for i, item in enumerate(mass_raw) if item not in set(mass_de)]
    if(len(diff_index)>1):
        # print("i am in wrong if")
        intensity_diff = list(itemgetter(*diff_index)(intensity_raw))

        return(intensity_diff)
    elif(len(diff_index)==1):
        # print("i am in right if")
        return([float(intensity_raw[diff_index[0]])])
    else:
        return(-1)



# print("i am spectra operation")

