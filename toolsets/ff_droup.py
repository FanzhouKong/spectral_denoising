import time

if __name__ == '__main__':
    freeze_support()

import pandas as pd
import numpy as np
import os
import toolsets.spectra_operations as so
from toolsets.search import string_search, quick_search_values, num_search
import toolsets.helpers as helpers
import pybaselines
import pymzml
from scipy.signal import find_peaks
# from multiprocessing import Pool, freeze_support
from matplotlib import rcParams
import matplotlib.pyplot as plt
from toolsets.raw_data_scaffold import read_mzml
import seaborn as sns
from toolsets.constants import iso_steps
iso_step_C = iso_steps['C']
iso_step_halo = iso_steps['Cl']
iso_step_H = iso_steps['H']

def feature_finding(mix, mzml_dir):
    '''
    wrapper function for get feature
    :return: data frame with all features
    '''
    ms1, ms2 = read_mzml(mix, parent_dir=mzml_dir)
    features_all = get_feature(ms1, ms2)
    return features_all
def get_feature(ms1_raw, ms2_raw, filter = True):
    ms1 = ms1_raw.copy()
    ms2 = ms2_raw.copy()
    max_length = len(ms2)
    pmz_list =np.zeros(max_length)
    rt_apex_list = np.zeros(max_length)
    rt_start_list = np.zeros(max_length)
    rt_end_list = np.zeros(max_length)
    rt_offset_list = np.zeros(max_length)
    mz_offset_list = np.zeros(max_length)
    reci_snr_list = np.zeros(max_length)
    peak_purity_list =np.zeros(max_length)
    msms_list = [None]*max_length
    peak_apex_intensity_list =np.zeros(max_length)
    peak_ms1_scan_range_list = [None]*max_length
    ms2_scan_idx_range_list = [None]*max_length
    ms2_scan_idx_list = np.zeros(max_length, dtype=int)
    mix_list = [None]*max_length
    current_idx = 0
    mass_sorted, intensity_sorted, index_sorted, rt_list = build_index(ms1)
    start = time.time()
    while(len(ms2))>0:
        ms2_pmz_sorted = ms2.sort_values(by = 'ms1_pmz', ascending=True)
        ms2_rt_sorted = ms2.sort_values(by = 'ms1_rt', ascending=True)
        seed_ms2 = ms2.iloc[np.argmax(ms2['ms1_precursor_intensity'])]
        pmz = seed_ms2['ms1_pmz']
        intensity_list = flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted)

        # intensity_list_smoothed = moving_average(intensity_list, n_neighbor=n_neighbor)
        baseline = pybaselines.smooth.snip(
            intensity_list, max_half_window=30, decreasing=True, smooth_half_window = 5
        )[0]
        baseline = [0 if x <0 else x for x in baseline]
        intensity_list_sub =[x-y for x, y in zip(intensity_list,baseline)]
        intensity_list_sub = [0 if x <0 else x for x in intensity_list_sub]
        n_neighbor = 2
        intensity_list_sub = moving_average(intensity_list_sub, n_neighbor=n_neighbor)
        peak_list = get_peaks(intensity_list_sub)
        target_peak_idx = find_most_close(peak_list, rt_list, seed_ms2['ms1_rt'], return_index=True)
        if len(peak_list) >0 and target_peak_idx != -1:
            target_peak = connect_peaks(peak_list, target_peak_idx, intensity_list_sub, rt_list)
            noise_level = np.median(baseline[target_peak[0]:target_peak[2]+1])
            rt_start, rt_end = rt_list[target_peak[0]], rt_list[target_peak[2]]
            rt_apex, int_apex = get_centroid(target_peak, rt_list, intensity_list_sub)
            if rt_apex != rt_apex:

                # return(seed_ms2, ms2)
                rt_apex = rt_list[target_peak[1]]
            reci_snr = noise_level/intensity_list[peak_list[target_peak_idx][1]]
            mapped_ms2_index=find_scan_fast(pmz, rt_start-5e-07, rt_end+5e-07, ms2_pmz_sorted, ms2_rt_sorted)
            mapped_ms2 = ms2.loc[mapped_ms2_index]
            feature_pmz = (mapped_ms2['ms1_pmz']*mapped_ms2['ms1_precursor_intensity']).sum()/(mapped_ms2['ms1_precursor_intensity'].sum())
            pmz_list[current_idx]=feature_pmz
            rt_apex_list[current_idx]=rt_apex
            rt_start_list[current_idx]=rt_start
            rt_end_list[current_idx]=rt_end
            rt_offset_list[current_idx]=(abs(rt_apex-seed_ms2['rt']))
            mz_offset_list[current_idx]=(abs(feature_pmz-seed_ms2['precursor_mz']))
            reci_snr_list[current_idx]=(reci_snr)
            peak_purity_list[current_idx]=(seed_ms2['peak_purity'])
            msms_list[current_idx]=(seed_ms2['peaks'])
            peak_apex_intensity_list[current_idx]=(int_apex)
            peak_ms1_scan_range_list[current_idx]=(list(target_peak))
            ms2_scan_idx_range_list[current_idx]=(list(mapped_ms2['scan_idx']))
            ms2_scan_idx_list[current_idx]=(seed_ms2['scan_idx'])
            mix_list[current_idx]=(seed_ms2['mix'])

            ms2.drop(mapped_ms2_index, inplace=True)
            current_idx = current_idx+1
        else:
            ms2 = ms2[ms2['scan_idx']!=seed_ms2['scan_idx']]

        # print(len(ms2))
        # break
    end = time.time()
    # print(end-start)
    features = pd.DataFrame(zip(pmz_list,
                                rt_apex_list,
                                peak_apex_intensity_list,
                                rt_start_list,
                                rt_end_list,
                                rt_offset_list,
                                mz_offset_list,
                                reci_snr_list,
                                peak_purity_list,
                                msms_list,
                                peak_ms1_scan_range_list,
                                ms2_scan_idx_range_list,
                                ms2_scan_idx_list,
                                mix_list),
                            columns=['precursor_mz', 'rt_apex', 'ms1_intensity',
                                     'rt_start', 'rt_end','rt_offset', 'mz_offset', 'reci_snr','peak_purity',
                                     'msms', 'ms1_scan_range', 'ms2_scan_idx_range_list',
                                     'ms2_scan_idx','mix']
                            )
    # features.sort_values(by = 'precursor_mz')
    features = features.iloc[0:current_idx]
    if filter == True:
        features = features[features['reci_snr']<1/3]
    iso_state= []
    is_dut = []
    is_halo = []
    for index, row in features.iterrows():
        iso_states_temp = _determine_iso_state(row, mass_sorted, intensity_sorted, index_sorted, )
        iso_state.append(iso_states_temp[0])
        is_dut.append(iso_states_temp[1])
        is_halo.append(iso_states_temp[2])
    features.insert(5, 'iso_state', iso_state)
    features.insert(6, 'is_dut', is_dut)
    features.insert(7, 'is_halo', is_halo)
    peaks_string = [None]*len(features)
    cur_idx = 0
    for index, row in features.iterrows():
        mass = row['msms'].T[0]
        intensity = row['msms'].T[1]
        peaks_string[cur_idx]=so.pack_spectra(mass, intensity)
        cur_idx = cur_idx+1
    features['msms']=peaks_string
    return(features)
def connect_peaks(peak_list, target_peak_idx, intensity_list, rt_list, split_tolerance = 1.3):
    apex_peak = list(peak_list[target_peak_idx])
    for i in range(1, target_peak_idx+1):
        left_peak = peak_list[target_peak_idx-i]
        # break
        if apex_peak[0]-left_peak[2]<=2 and intensity_list[apex_peak[0]]>min(intensity_list[apex_peak[1]], intensity_list[left_peak[1]])/split_tolerance:
            apex_peak[0] = left_peak[0]
            apex_peak[1]=np.argmax(intensity_list[apex_peak[0]:apex_peak[2]])+apex_peak[0]
        else:
            break
    for i in range(1, len(peak_list)-target_peak_idx):
        right_peak = peak_list[target_peak_idx+i]
        # break
        if right_peak[0]-apex_peak[2]<=2 and intensity_list[apex_peak[2]]>min(intensity_list[apex_peak[1]], intensity_list[right_peak[1]])/split_tolerance:
            apex_peak[2] = right_peak[2]
            apex_peak[1]=np.argmax(intensity_list[apex_peak[0]:apex_peak[2]])+apex_peak[0]
        else:
            break
    return(tuple(apex_peak))
def auto_EIC(pmz, ms1, mass_tolerance = 0.005, show_eic = False):
    mass_sorted, intensity_sorted, index_sorted, rt_list = build_index(ms1)
    intensity_list = flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted, mass_error=mass_tolerance)
    if show_eic == True:
        EIC(rt_list, intensity_list)
        return()
    return rt_list, intensity_list
def EIC(rt_list, intensity_list,
        # parent_dir= None, if_mix = False,
        base_line_level = -1,
        base_line_series = None,
        rt_start = -1, rt_end = -1, adjusted_height = -1, vlines_location_1 = [], vlines_location_2 = [],annotation = None,
        savepath = None, show =True):
    # print('tttt')
    fig, ax = plt.subplots(
        figsize = (10, 6)
    )
    rcParams.update({'font.size':8})
    ax= sns.lineplot(x = rt_list, y = intensity_list, label = 'EIC')
    if base_line_series is not None:
        if len(base_line_series)> len(rt_list):
            start_idx = int(np.floor(len(base_line_series)/2)-np.floor(len(rt_list)/2))
            end_idx = int(np.floor(len(base_line_series)/2)+np.floor(len(rt_list)/2)+1)
            res = base_line_series[start_idx:end_idx]
            ax = sns.lineplot(x = rt_list, y = res, color = 'orange', label = 'baseline')
        elif len(base_line_series)<len(rt_list):
            offset = int((len(rt_list)-len(base_line_series))/2)
            print(offset)
            rt_list = rt_list[offset:-offset]
            intensity_list = intensity_list[offset:-offset]
            print(len(base_line_series))
            print(len(rt_list))
            ax = sns.lineplot(x = rt_list, y = base_line_series, color = 'orange', label = 'baseline')
        elif len(base_line_series)==len(rt_list):
            ax = sns.lineplot(x = rt_list, y = base_line_series, color = 'orange', label = 'baseline')
    ax.set_ylim(0, np.max(intensity_list)+100)
    if rt_start != -1 and rt_end != -1:
        index_start = np.searchsorted(rt_list, rt_start,side = 'left')
        index_end = np.searchsorted(rt_list, rt_end,side = 'right')
        # rt_list = rt_list[index_start:index_end+1]

        ax.set_xlim(rt_start, rt_end)
        ax.set_ylim(0, adjusted_height)
        ax.set_ylim(0, np.max(intensity_list[index_start:index_end])*1.1)
    # else:
    #     ax.set_xlim(0, np.max(rt_list))
    if base_line_level != -1:
        plt.axhline(base_line_level, color = 'orange')
    if adjusted_height!= -1:
        ax.set_ylim(0, adjusted_height)
    if len(vlines_location_1)>0:
        for position in vlines_location_1:
            plt.axvline(x = position, color = 'red')
    if len(vlines_location_2)>0:
        for position in vlines_location_2:
            plt.axvline(x = position, color = 'green')
    ax.grid(False)
    ax.set_facecolor("none")
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if annotation is not None:
        ax.set_title(str(annotation))
    # ax.set(xticklabels=['Retention time (min)'], yticklabels = ['intensity'])
    fig.tight_layout()
    # ax.patch.set_linewidth(1)
    # ax.set_edgecolor("black")
    if savepath != None:
        plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'black')
    if show != True:
        plt.close()

from toolsets.search import quick_search_sorted
def find_scan_fast(pmz,rt_start, rt_end, ms2_pmz_sorted,ms2_rt_sorted,mass_tolerance = 0.005):
    # feature_mz_sorted = feature_table.sort_values(by = mz_column, ascending=True)
    # feature_rt_sorted = feature_table.sort_values(by = rt_column, ascending=True)
    scan_pmz_matched = quick_search_sorted(ms2_pmz_sorted, 'ms1_pmz', pmz-mass_tolerance, pmz+mass_tolerance)
    scan_rt_matched = quick_search_sorted(ms2_rt_sorted, 'ms1_rt', rt_start, rt_end)
    lst = list(set(scan_pmz_matched.index) & set(scan_rt_matched.index))
    return(lst)

def moving_average( intensity_list, n_neighbor = 3):
    # n_neighbor = 20
    bl  = []
    bl_rt = []
    for i in range(len(intensity_list)):
        if i < n_neighbor or len(intensity_list)-i<=n_neighbor:
            bl.append(intensity_list[i])
        else:
            neighbors = intensity_list[i-n_neighbor:i+n_neighbor]
            bl.append(np.mean(neighbors))
    # print(n_neighbor)

    return( bl)

from operator import itemgetter

from scipy.stats import pearsonr
def _determine_iso_state(row, mass_sorted, intensity_sorted, index_sorted, iso_range = 2, iso_corrs_threshold = 0.75,if_debug = False):

    ion_trace_temp = flash_eic(row['precursor_mz'], mass_sorted, intensity_sorted, index_sorted)
    index_start = row['ms1_scan_range'][0]
    index_end = row['ms1_scan_range'][2]+1
    index_apex = row['ms1_scan_range'][1]
    halo_check = _check_halogen(row['precursor_mz'], mass_sorted, intensity_sorted, index_sorted, index_start, index_apex, index_end, iso_corrs_threshold = iso_corrs_threshold)
    if_dut = False
    if_halo = False
    iso_state = 0
    if halo_check[0]==True:
        if_halo = halo_check[0]
        iso_state = halo_check[1]
        mono_isotopic_mass = row['precursor_mz']-iso_state/2*iso_step_halo
        ion_trace_mono = flash_eic(mono_isotopic_mass, mass_sorted, intensity_sorted, index_sorted)
        ion_trace_dut = flash_eic(mono_isotopic_mass-iso_step_H, mass_sorted, intensity_sorted, index_sorted)
        corr, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_dut[index_start:index_end])
        if corr == corr and corr >iso_corrs_threshold and ion_trace_dut[index_apex]<ion_trace_mono[index_apex]:
            if_dut = True
        return(iso_state, if_dut, if_halo)

    corrs = np.zeros(2*iso_range+1)
    iso_intensity = np.zeros(2*iso_range+1)
    for step in range(-iso_range,iso_range+1):
        ion_trace_iso = flash_eic(row['precursor_mz']+(step*iso_step_C), mass_sorted, intensity_sorted, index_sorted)
        # print(row['precursor_mz']+(step*iso_step))
        lst_idx = step+iso_range
        corr, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_iso[index_start:index_end])
        iso_intensity[lst_idx]=(ion_trace_iso[row['ms1_scan_range'][1]])
        if corr != corr:
            corrs[lst_idx]=(0)
        else:
            corrs[lst_idx]=(corr)
    if np.argmax(iso_intensity) == 0:
        if_dut = False
    elif corrs[np.argmax(iso_intensity)-1]>iso_corrs_threshold:
        if_dut= True
    else:
        if_dut= False
    # return(iso_range, corrs, iso_intensity)
    iso_state_temp = iso_range-np.argmax(iso_intensity)
    # print(np.argmax(iso_intensity))
    # print(len(iso_intensity))
    if np.argmax(iso_intensity)+1 == len(iso_intensity) and np.max(corrs)>iso_corrs_threshold:
        iso_state = iso_state_temp
    elif np.argmax(iso_intensity)+1 == len(iso_intensity):
        iso_state = np.NAN

    elif corrs[np.argmax(iso_intensity)+1]>iso_corrs_threshold:
        iso_state = iso_state_temp
        # print(corrs[np.argmax(iso_intensity)-1],iso_corrs_threshold)



    if if_debug == False:
        return(iso_state, if_dut, if_halo)
    else:
        return(corrs, iso_intensity)
def _check_halogen(pmz, mass_sorted, intensity_sorted, index_sorted, index_start, index_apex, index_end, iso_corrs_threshold = 0.7):
    iso_step = iso_step_halo
    is_halo = False
    iso_state = 0
    ion_trace_temp = flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted)
    ion_trace_C_left = flash_eic(pmz-iso_step_C, mass_sorted, intensity_sorted, index_sorted)
    ion_trace_C_right = flash_eic(pmz+iso_step_C, mass_sorted, intensity_sorted, index_sorted)
    ion_trace_halo_left = flash_eic(pmz-iso_step, mass_sorted, intensity_sorted, index_sorted)
    ion_trace_halo_right = flash_eic(pmz+iso_step, mass_sorted, intensity_sorted, index_sorted)

    corr_left, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_halo_left[index_start:index_end])
    corr_right, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_halo_right[index_start:index_end])

    if corr_left == corr_left and  corr_left>iso_corrs_threshold and ion_trace_C_left[index_apex] <ion_trace_temp[index_apex] and ion_trace_C_left[index_apex] <ion_trace_halo_left[index_apex]:
        is_halo = True
        search_left = True
    elif corr_right == corr_right and corr_right>iso_corrs_threshold and ion_trace_C_right[index_apex] <ion_trace_temp[index_apex] and ion_trace_C_right[index_apex] <ion_trace_halo_right[index_apex]:
        is_halo = True
        search_left = False
    else:
        return(is_halo, iso_state)
    if search_left == False:
        return(is_halo, 0)
    else:
        round = 1
        while(search_left==True):
            pmz = pmz-iso_step
            ion_trace_temp = flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted)
            ion_trace_C_left = flash_eic(pmz-iso_step_C, mass_sorted, intensity_sorted, index_sorted)
            ion_trace_halo_left = flash_eic(pmz-iso_step, mass_sorted, intensity_sorted, index_sorted)
            corr_left, _ = pearsonr(ion_trace_temp[index_start:index_end], ion_trace_halo_left[index_start:index_end])
            if corr_left == corr_left and corr_left>iso_corrs_threshold and ion_trace_C_left[index_apex] <ion_trace_temp[index_apex] and ion_trace_C_left[index_apex] <ion_trace_halo_left[index_apex]:
                round = round+1
            else:
                search_left = False
        return(is_halo, round*2)
import itertools
def build_index(ms1):
    # ms1.reset_index(inplace = True, drop = True)
    mass_nested = [None]*len(ms1)
    intensity_nested = [None]*len(ms1)
    index_nested = [None]*len(ms1)
    rt_list = np.zeros(len(ms1))
    for index, row in (ms1.iterrows()):
        mass_temp, intensity_temp = row['peaks'].T
        mass_nested[index]=(mass_temp)
        intensity_nested[index]=(intensity_temp)
        index_nested[index]=([index]*len(mass_temp))
        rt_list[index]=(row['rt'])
    mass_flatten = np.array(list(itertools.chain.from_iterable(mass_nested)))
    intensity_flatten = np.array(list(itertools.chain.from_iterable(intensity_nested)))
    index_flatten = np.array(list(itertools.chain.from_iterable(index_nested)))
    order = np.argsort(mass_flatten)
    mass_sorted = mass_flatten[order]
    intensity_sorted = intensity_flatten[order]
    index_sorted = index_flatten[order]
    return(mass_sorted, intensity_sorted, index_sorted, rt_list)
# def build_index(ms1):
#     '''
#     deprecated
#     '''
#     mass_flatten = []
#     intensity_flatten = []
#     index_flatten = []
#     rt_list = []
#     for index, row in (ms1.iterrows()):
#         # mass_temp, intensity_temp = so.break_spectra(row['peaks'])
#         mass_temp, intensity_temp = row['peaks'].T
#         mass_flatten.extend(mass_temp)
#         intensity_flatten.extend(intensity_temp)
#         index_flatten.extend([index]*len(mass_temp))
#         rt_list.append(row['rt'])
#     mass_sorted, intensity_sorted, index_sorted = zip(*sorted(zip(mass_flatten, intensity_flatten, index_flatten)))
#     mass_sorted = np.array(mass_sorted)
#     intensity_sorted= np.array(intensity_sorted)
#     index_sorted = np.array(index_sorted)
#     return(mass_sorted, intensity_sorted, index_sorted, rt_list)

def flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted, mass_error=0.005):
    # print('ttt')
    # index_start, index_end = mass_sorted.searchsorted([pmz-mass_error, pmz+mass_error+1E-9])
    index_start = np.searchsorted(mass_sorted, pmz-mass_error, side = 'left')
    index_end = np.searchsorted(mass_sorted, pmz+mass_error, side = 'right')
    index_range = index_sorted[index_start:index_end]

    intensity_range = intensity_sorted[index_start:index_end]

    intensity_list = np.zeros(np.max(index_sorted)+1)
    for idx in range(0,len(index_range)):
        intensity_list[index_range[idx]]= intensity_list[index_range[idx]]+intensity_range[idx]
    return(intensity_list)


def find_most_close(peaks, rt_list, rt_apex, return_index = False):
    if len(peaks)==0:
        return(-1)
    offset = []
    # peak_contain = []
    for peak in peaks:
        # offset.append(abs(rt_apex - rt_list[peak[1]]))
        if rt_apex >=rt_list[peak[0]] and rt_apex <= rt_list[peak[2]]:
            offset.append(abs(rt_apex - rt_list[peak[1]]))
        else:
            offset.append(np.NAN)
    # return(offset)
    try:
        return(np.nanargmin(offset))
    except:
        return(-1)




def get_peaks(intensity_list: np.ndarray) -> list:
    """Detects peaks in an array.

    Args:
        int_array (np.ndarray): An array with intensity values.

    Returns:
        list: A regular Python list with all peaks.
            A peak is a triplet of the form (start, center, end)

    """
    apex, _ = find_peaks(intensity_list)
    peak_list = []
    for cur_apex_idx in apex:
        peak_list.append(get_edges(intensity_list, cur_apex_idx))
    return(peak_list)

def get_edges(intensity_list, cur_apex_idx):
    intensity_list = np.array(intensity_list)
    gradient = np.diff(intensity_list)
    left_edge_idx = cur_apex_idx
    right_edge_idx = cur_apex_idx
    while left_edge_idx>0:
        if intensity_list[left_edge_idx-1]<=intensity_list[left_edge_idx] and intensity_list[left_edge_idx]>0:
            left_edge_idx = left_edge_idx-1
        else:
            break
    while right_edge_idx <len(intensity_list)-1:
        if intensity_list[right_edge_idx+1]<=intensity_list[right_edge_idx] and intensity_list[right_edge_idx]>0:
            right_edge_idx = right_edge_idx+1
        else:
            break
    # while left_edge_idx>0:
    #     if gradient[left_edge_idx-1]>0 and intensity_list[left_edge_idx]>0:
    #         left_edge_idx = left_edge_idx-1
    #     else:
    #         break
    # while right_edge_idx < len(intensity_list)-1:
    #     if gradient[right_edge_idx+1]<0 and intensity_list[left_edge_idx]>0:
    #         right_edge_idx = right_edge_idx+1
    #     else:
    #         break
    return((left_edge_idx, cur_apex_idx, right_edge_idx))



def get_centroid(
        peak: tuple,
        mz_array: np.ndarray,
        int_array: np.ndarray
) -> tuple:
    """Wrapper to estimate centroid center positions.

    Args:
        peak (tuple): A triplet of the form (start, center, end)
        mz_array (np.ndarray): An array with mz values.
        int_array (np.ndarray): An array with intensity values.

    Returns:
        tuple: A tuple of the form (center, intensity)
    """
    start, center, end = peak
    mz_int = np.sum(int_array[start : end+1])
    mz_apex = int_array[center]

    peak_size = end - start - 1

    if peak_size == 1:
        mz_cent = mz_array[center]
    elif peak_size == 2:
        mz_cent = (
                          mz_array[start + 1] * int_array[start + 1]
                          + mz_array[end - 1] * int_array[end - 1]
                  ) / (int_array[start + 1] + int_array[end - 1])
        # print('peak_size = 2')
    else:
        mz_cent = gaussian_estimator(peak, mz_array, int_array)

    # return mz_cent, mz_int
    # return mz_cent, mz_int
    return mz_cent, mz_apex
def gaussian_estimator(
        peak: tuple,
        mz_array: np.ndarray,
        int_array: np.ndarray
) -> float:
    """Three-point gaussian estimator.

    Args:
        peak (tuple): A triplet of the form (start, center, end)
        mz_array (np.ndarray): An array with mz values.
        int_array (np.ndarray): An array with intensity values.

    Returns:
        float: The gaussian estimate of the center.
    """
    start, center, end = peak

    m1, m2, m3 = mz_array[center - 1], mz_array[center], mz_array[center + 1]
    i1, i2, i3 = int_array[center - 1], int_array[center], int_array[center + 1]
    # print(m1,m2,m3)
    if i1 == 0:  # Case of sharp flanks
        m = (m2 * i2 + m3 * i3) / (i2 + i3)
    elif i3 == 0:
        m = (m1 * i1 + m2 * i2) / (i1 + i2)
    else:
        l1, l2, l3 = np.log(i1), np.log(i2), np.log(i3)
        m = (
                ((l2 - l3) * (m1 ** 2) + (l3 - l1) * (m2 ** 2) + (l1 - l2) * (m3 ** 2))
                / ((l2 - l3) * (m1) + (l3 - l1) * (m2) + (l1 - l2) * (m3))
                * 1
                / 2
        )

    return m
def centroid_data(
        mz_array: np.ndarray,
        int_array: np.ndarray
) -> tuple:
    """Estimate centroids and intensities from profile data.

    Args:
        mz_array (np.ndarray): An array with mz values.
        int_array (np.ndarray): An array with intensity values.

    Returns:
        tuple: A tuple of the form (mz_array_centroided, int_array_centroided)
    """
    peaks = get_peaks(int_array)

    mz_array_centroided = np.zeros(len(peaks))
    int_array_centroided = np.zeros(len(peaks))


    for i in range(len(peaks)):
        mz_, int_ = get_centroid(peaks[i], mz_array, int_array)
        mz_array_centroided[i] = mz_
        int_array_centroided[i] = int_

    return mz_array_centroided, int_array_centroided
import plotly.express as px
#
# def feature_finding(ms1, ms2):
#     features = pd.DataFrame()
#     if len(ms2)==0:
#         return(features)
#     # peak_purity_tolerance = ms2['peak_purity'].describe()['mean']-ms2['peak_purity'].describe()['std']*2
#     # ms2 = quick_search_values(ms2, 'peak_purity', peak_purity_tolerance, 1)
#     bins = get_mz_bin(ms2)
#     with Pool(processes=6) as pool:
#         results= pool.starmap(get_feature, zip(repeat(ms2), bins, repeat(ms1)))
#     # print('i exited pool okay')
#
#     for result in results:
#         features = pd.concat([features, result], ignore_index=True)
#     if len(features)==0:
#         return features
#     features.sort_values(by = 'rt_offset', ascending=True, inplace=True)
#     features.drop_duplicates(subset=['ms2_range_idx'], keep = 'first',inplace=True, )
#     features.reset_index(inplace=True, drop=True)
#     return(features)
#     # for feature

    # pass

# def get_feature(ms2,current_pmz_bin,ms1):
#     # print('tttt')
#     current_cluster = quick_search_values(ms2, 'ms1_pmz', value_start=current_pmz_bin-ms1_tolerance, value_end=current_pmz_bin+ms1_tolerance)
#     # print('i am in new!')
#     features_all = pd.DataFrame()
#
#     rt_list, intensity_list = get_EIC_list(ms1, current_pmz_bin, step=ms1_tolerance)
#     n_neighbor = 2
#     intensity_list = moving_average(intensity_list, n_neighbor=n_neighbor)
#     baseline = pybaselines.smooth.snip(
#         intensity_list, max_half_window=30, decreasing=True, smooth_half_window = 5
#     )[0]
#     baseline = [0 if x <0 else x for x in baseline]
#     intensity_list_raw = intensity_list.copy()
#     intensity_list = [x-y for x, y in zip(intensity_list_raw, baseline)]
#     intensity_list=[0 if x <0 else x for x in intensity_list]
#     peak_list = get_peaks(intensity_list)
#
#     while len(current_cluster)>0:
#
#         seed_ms2 = current_cluster.iloc[np.argmax(current_cluster['ms1_precursor_intensity'])]
#
#         target_peak_idx = find_most_close(peak_list, rt_list, seed_ms2['ms1_rt'], return_index=True)
#         if len(peak_list) >0 and target_peak_idx != -1:
#
#             target_peak = connect_peaks(peak_list, target_peak_idx, intensity_list, rt_list)
#             rt_start, rt_end = rt_list[target_peak[0]], rt_list[target_peak[2]]
#             rt_apex, int_apex = get_centroid(target_peak, rt_list, intensity_list)
#             if rt_apex != rt_apex:
#                 rt_apex = rt_list[target_peak[1]]
#             baseline_median = np.median(baseline[target_peak[0]:target_peak[2]+1])
#             if baseline_median == 0:
#                 baseline_median =np.mean(baseline)
#             snr = (intensity_list_raw[target_peak[1]]/baseline_median)
#
#             # current_cluster_rt_sorted = current_cluster.sort_values(by = 'ms1_rt', ascending=True)
#             ms2_under_peaks  = quick_search_values(current_cluster, 'ms1_rt', value_start=rt_start, value_end=rt_end)
#             if len(ms2_under_peaks)>0 and snr >3:
#                 ms2_under_peaks['rt']=rt_apex
#                 ms2_under_peaks['rt_offset']=abs(rt_apex-seed_ms2['rt'])
#                 ms2_under_peaks['rt_start']=rt_start
#                 ms2_under_peaks['rt_end']=rt_end
#                 ms2_under_peaks['precursor_mz']=(ms2_under_peaks['ms1_pmz']*ms2_under_peaks['ms1_precursor_intensity']).sum()/(ms2_under_peaks['ms1_precursor_intensity'].sum())
#
#                 ms2_under_peaks['peak_apex_intensity']=int_apex
#                 ms2_under_peaks['peaks_all']=ms2_under_peaks.apply(lambda x: list(ms2_under_peaks['peaks']), axis=1)
#                 ms2_under_peaks['ms1_intensities']=ms2_under_peaks.apply(lambda x: list(ms2_under_peaks['ms1_precursor_intensity']), axis=1)
#                 key_peak_idx = np.argmax(ms2_under_peaks['ms1_precursor_intensity'])
#                 key_peak = ms2_under_peaks.iloc[key_peak_idx]['peaks']
#                 ms2_under_peaks['msms']=key_peak
#                 ms2_under_peaks['wa_msms']=so.weighted_average_spectra(ms2_under_peaks, weight_col= 'ms1_precursor_intensity')
#                 ms2_under_peaks['snr']=snr
#                 ms2_under_peaks['peak_range_idx'] = ms2_under_peaks.apply(lambda x: list(target_peak), axis=1)
#                 ms2_under_peaks['ms2_range_idx'] = ms2_under_peaks.apply(lambda x: list(ms2_under_peaks['scan_idx'].unique()), axis=1)
#                 ms2_under_peaks['pmz_bin']=current_pmz_bin
#                 ms2_under_peaks.drop(['scan_idx', 'cycle', 'isolation_window','peaks', 'ms1_pmz', 'ms1_precursor_intensity', 'peak_purity', 'mz_offset', 'ms1_rt','ms1_intensities','peaks_all'], axis =1, inplace = True)
#                 current_cluster.drop(ms2_under_peaks.index, inplace = True)
#                 features_all=pd.concat([features_all,pd.DataFrame([ms2_under_peaks.iloc[0]]) ], axis=0)
#             else:
#                 current_cluster = string_search(current_cluster, 'scan_idx', seed_ms2.scan_idx, reverse=True)# skip this seed ms2
#         else:
#             current_cluster = string_search(current_cluster, 'scan_idx', seed_ms2.scan_idx, reverse=True)# skip this seed ms2
#            # break
#     if len(features_all)>0:
#         features_all.sort_values(by = 'rt_offset', ascending=True, inplace=True)
#         features_all.drop_duplicates(subset='rt', keep='first', inplace=True)
#     return (features_all)
    # if len(features_all)>0:
    #     features_tidy = pd.DataFrame()
    #     for rt in features_all['rt'].unique():
    #         feature_rt = string_search(features_all, 'rt', rt)
    #         feature_rt.sort_values(by = 'rt_offset', inplace= True, ignore_index=True, ascending=True)
    #         # features_tidy= features_tidy.append(feature_rt.iloc[0])
    #         features_tidy = pd.concat([features_tidy,pd.DataFrame([feature_rt.iloc[0]])], axis =0)
    #     features_tidy.reset_index(inplace=True, drop=True)
    #
    #     return(features_tidy)
    # else:
    #     return(features_all)




# def get_EIC_list(ms1, pmz, step = 0.005):
#     rt_list = []
#     intensity_list = []
#
#     for index, row in ms1.iterrows():
#         rt_list.append(row['rt'])
#         intensity_list.append(_extract_ms1_intensity(row['peaks'], mz_lower=pmz-step, mz_upper=pmz+step))
#     return(rt_list, intensity_list)
# def auto_EIC(mix, parent_dir,pmz, vlines_location_1=[], vlines_location_2=[] , rt_start = -1, rt_end = -1, rt_max = 12):
#     ms1, ms2 = process_mzml(mix, parent_dir, if_mix=True, with_ms1=True, rt_max = rt_max)
#     rt_list, intensity_list = get_EIC_list(ms1, pmz)
#     EIC(rt_list, intensity_list, vlines_location_1=vlines_location_1,vlines_location_2 = vlines_location_2, rt_start=rt_start, rt_end=rt_end)

    # peaklist = []
    # gradient = np.diff(int_array)
    # start, center, end = -1, -1, -1

    # for i in range(len(gradient)):

    #     grad = gradient[i]

    #     if (end == -1) & (center == -1):  # No end and no center yet
    #         if grad <= 0:  # If decreasing, move start point
    #             start = i
    #         else:  # If increasing set as center point
    #             center = i

    #     if (end == -1) & (
    #         center != -1
    #     ):  # If we have a centerpoint and it is still increasing set as new center
    #         if grad >= 0:
    #             center = i
    #         else:  # If it is decreasing set as endpoint
    #             end = i

    #     if end != -1:  # If we have and endpoint and it is going down
    #         if grad < 0:
    #             end = i  # Set new Endpoint
    #         else:  # if it stays the same or goes up set a new peak
    #             peaklist.append((start + 1, center + 1, end + 1))
    #             start, center, end = end, -1, -1  # Reset start, center, end

    # if end != -1:
    #     peaklist.append((start + 1, center + 1, end + 1))

    # return peaklist

# def get_mz_bin(ms2, step = 0.005):
#     mass_working = np.array(ms2['ms1_pmz'])
#     intensity_working = np.array(ms2['ms1_precursor_intensity'])
#     pmz_bin = []
#     while len(mass_working)>0:
#         idx = np.argmax(intensity_working)
#         start_index, end_index = mass_working.searchsorted([mass_working[idx]-step, mass_working[idx]+step])
#         pmz_bin.append(mass_working[idx])
#         mass_working = np.delete(mass_working, np.arange(start_index, end_index))
#         intensity_working = np.delete(intensity_working, np.arange(start_index, end_index))
#
#     return(pmz_bin)
    # valley, _ = find_peaks(intensity_list_neg)
    # with contextlib.redirect_stdout(None):
    #    fp_model = fp.findpeaks(lookahead=1,interpolate=10)
    #    # fit
    #    results=fp_model.fit(intensity_temp)['df']
    # results_apex = string_search(results, 'peak', True, reset_index=False)
    # left_edge = []
    # right_edge = []
    # apex = []
    # for index, row in results_apex.iterrows():
    #    apex.append(row['x'])
    #    for i in range(1,row['x']+1):
    #       if results.iloc[row['x']-i]['valley']==True:
    #          left_edge_temp = results.iloc[row['x']-i]['x']
    #          break
    #    left_edge.append((left_edge_temp))
    #    for i in range(1, len(results)-row['x']):
    #       if results.iloc[row['x']+i]['valley']==True:
    #          right_edge_temp =results.iloc[row['x']+i]['x']
    #          break
    #    right_edge.append(right_edge_temp)
    # peak_list = []
    #
    # for idx in range(0, len(apex)):
    #
    #     peak_list.append((left_edge[idx], apex[idx], right_edge[idx]))
    # return(peak_list)
# def get_istd_info_all(istd_info, file_list, mzml_dir):
#
#     with Pool(processes=6) as pool:
#         results= pool.starmap(get_istd_info, zip(repeat(istd_info), file_list, repeat(mzml_dir)))
#     # print('i exited pool okay')
#     # features = pd.DataFrame()
#     # for result in results:
#     #     features = pd.concat([features, result], ignore_index=True)
#     # features.sort_values(by = 'rt_offset', ascending=True, inplace=True)
#     # features.drop_duplicates(subset=['ms2_range_idx'], keep = 'first',inplace=True, )
#     # features.reset_index(inplace=True, drop=True)
#     return(results)
# def get_istd_info(istd_info, filename, mzml_dir):
#
#     ms1, ms2 = process_mzml(mzml_path=os.path.join(mzml_dir, filename+'.mzML'), rt_max=5)
#     intensity_file = []
#     rt_file = []
#     for index, row in istd_info.iterrows():
#         rt_temp, intensity_temp = get_ms1_feature(row['pmz'],row['rt_suggested'], ms1)
#         rt_file.append(rt_temp)
#         intensity_file.append(intensity_temp)
#     return(filename, rt_file, intensity_file)
# def get_ms1_feature(pmz,  rt,ms1):
#     rt_list, intensity_list = get_EIC_list(ms1, pmz, step = ms1_tolerance)
#     n_neighbor = 2
#     intensity_list = moving_average(intensity_list, n_neighbor=n_neighbor)
#     rt_list = rt_list[n_neighbor:-n_neighbor]
#     baseline = pybaselines.smooth.snip(
#         intensity_list, max_half_window=30, decreasing=True, smooth_half_window = 5
#     )[0]
#     baseline = [0 if x <0 else x for x in baseline]
#     intensity_list_raw = intensity_list.copy()
#     intensity_list = [x-y for x, y in zip(intensity_list_raw, baseline)]
#     intensity_list=[0 if x <0 else x for x in intensity_list]
#     peak_list = get_peaks(intensity_list)
#     target_peak_idx = find_most_close(peak_list, rt_list,rt, return_index=True)
#     candidate_intensity = []
#     #
#     for peak in peak_list:
#         # offset.append(abs(rt_apex - rt_list[peak[1]]))
#         if rt_list[peak[1]]>rt-20/60 and rt_list[peak[1]]<rt+20/60:
#             # candidate_peak.append(peak)
#             # candidate_idx.append(peak[1])
#             candidate_intensity.append(intensity_list[peak[1]])
#         else:
#             candidate_intensity.append(np.NAN)
#     # print(apex_idx)return(apex_idx)
#     if np.isnan(candidate_intensity).all()==False:
#         apex_idx = np.nanargmax(candidate_intensity)
#
#         target_peak_idx = peak_list[apex_idx]
#
#
#         rt_apex, int_apex = get_centroid([target_peak_idx[1]-1,target_peak_idx[1], target_peak_idx[1]+1], rt_list, intensity_list)
#     else:
#         rt_apex = np.NAN
#         int_apex = np.NAN
#     return(rt_apex, int_apex)