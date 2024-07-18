import toolsets.raw_data_scaffold as rds
import numpy as np
import pandas as pd
from toolsets.search import string_search, quick_search_sorted, quick_search_values
import toolsets.helpers as helpers
import pybaselines
from tqdm import tqdm
import os
from scipy.signal import find_peaks
import toolsets.helpers as helpers
import toolsets.spectra_operations as so
from collections import Counter
def flash_bpc(ms1):
    bpc_list = np.zeros(len(ms1))
    idx = 0
    for index, row in (ms1.iterrows()):
        mass_temp, intensity_temp = row['peaks'].T
        bpc_list[idx] = bpc_list[idx]+np.max(intensity_temp)
        idx = idx+1
    return bpc_list
# def find_boundry(mass_sorted, intensity_sorted, tolerance = 30000):
#     mass_sorted_indexing = mass_sorted.copy()
#     intensity_sorted_indexing = intensity_sorted.copy()
#     left_boundry = []
#     right_boundry = []
#     check_len = len(intensity_sorted_indexing)
#     while np.max(intensity_sorted_indexing)>tolerance:
#         seed_idx = np.argmax(intensity_sorted_indexing)
#         #search left
#         # for l in np.flip(np.arange(1, seed_idx)):
#         #     diff = mass_sorted_indexing[l]-mass_sorted_indexing[l-1]
#         #     if diff >0.005:
#         #         break
#         # for r in np.arange(seed_idx, len(intensity_sorted_indexing)-1):
#         #     diff = mass_sorted_indexing[r+1] - mass_sorted_indexing[r]
#         #     if diff >0.005:
#         #         break
#         # if l != r:
#         #     left_boundry.append(mass_sorted_indexing[l])
#         #     right_boundry.append(mass_sorted_indexing[r])
#         indices_to_drop = np.arange(l, r+1)
#         mask = np.ones(len(mass_sorted_indexing), dtype=bool)
#         mask[indices_to_drop] = False
#         mass_sorted_indexing = mass_sorted_indexing[mask]
#         intensity_sorted_indexing = intensity_sorted_indexing[mask]
#         if check_len == len(intensity_sorted_indexing):
#             print('not shrinking')
#             return(mass_sorted[seed_idx])
#         else:
#             check_len = len(intensity_sorted_indexing)
#
#     return left_boundry, right_boundry
# def find_nearest_exceeding_indices(arr, index, number):
#     left_index = None
#     right_index = None
#
#     # Search to the left
#     for i in range(index, -1, -1):
#         if arr[i] > number:
#             left_index = i
#             break  # Stop as soon as you find the first exceeding number on the left
#
#     # Search to the right
#     for i in range(index, len(arr)):
#         if arr[i] > number:
#             right_index = i
#             break  # Stop as soon as you find the first exceeding number on the right
#
#     return left_index, right_index
def find_feature(mzml_name, mzml_dir, intensity_threshold = 30000, n_neighbor = 2):
    ms1, ms2 = rds.read_mzml(mzml_name, mzml_dir)
    mass_sorted, intensity_sorted, index_sorted, rt_list = build_index(ms1)
    feature_temp = get_features(mass_sorted, intensity_sorted, index_sorted, rt_list, intensity_threshold= intensity_threshold, n_neighbor=n_neighbor)
    feature_temp = deduplication(feature_temp)
    feature_temp_ms2 = map_ms2(feature_temp, ms2)
    return(feature_temp_ms2)
def deduplication(ms1_features):
    ms1_features['pmz_offset']=abs(ms1_features['precursor_mz']-ms1_features['eic_center'])
    ms1_features['anchor_scan']=ms1_features['ms1_scan_range'].apply(lambda x: x[1])
    masks = ms1_features.groupby(['anchor_pmz', 'anchor_scan']).pmz_offset.transform(min)
    ms1_deduplicated = ms1_features.loc[ms1_features.pmz_offset == masks]
    ms1_deduplicated.drop(columns = ['anchor_scan', 'anchor_pmz'], inplace = True)
    return ms1_deduplicated
def get_seed_mass(seed_mass, mass_sorted):
    if seed_mass<400:
        pre_mass_error = 0.004
    else:
        pre_mass_error = seed_mass * 10 / 1E6
    left_idx, right_idx = mass_sorted.searchsorted([seed_mass-pre_mass_error, seed_mass+pre_mass_error])
    mass_error = np.std(mass_sorted[left_idx:right_idx])*3
    if mass_error <pre_mass_error:
        mass_error = pre_mass_error
    left_idx, right_idx = mass_sorted.searchsorted([seed_mass-mass_error, seed_mass+mass_error])
    seed_mass = np.median(mass_sorted[left_idx:right_idx])
    left_idx, right_idx = mass_sorted.searchsorted([seed_mass-mass_error, seed_mass+mass_error])
    return(seed_mass, mass_error, left_idx, right_idx)
def is_in_ranges(number, centers, scopes):
    """
    Check if a given number falls within any of the ranges defined by centers and scopes.

    :param number: The number to check.
    :param centers: Array of range centers.
    :param scopes: Array of range scopes (half-widths).
    :return: True if the number is in any range, False otherwise.
    """
    # Calculate start and end points of each range
    starts = centers - scopes
    ends = centers + scopes

    # Check if the number is in any range
    in_range = np.any((number >= starts) & (number <= ends))

    return in_range
def get_features(mass_sorted, intensity_sorted, index_sorted, rt_list, intensity_threshold = 1000, n_neighbor = 2):
    pmz = np.zeros(len(mass_sorted))

    rt = np.zeros(len(mass_sorted))
    rt_start =np.zeros(len(mass_sorted))
    rt_end =np.zeros(len(mass_sorted))
    eic_center = np.zeros(len(mass_sorted))
    anchor_pmzs = np.zeros(len(mass_sorted))
    eic_offset = np.zeros(len(mass_sorted))
    apex_intensity = np.zeros(len(mass_sorted))
    n_scans = np.zeros(len(mass_sorted))
    peak_range = [None]*(len(mass_sorted))
    ms1_idx = [None]*(len(mass_sorted))
    reci_snrs = np.zeros((len(mass_sorted)))
    counter = 0

    order = np.argsort(intensity_sorted)[::-1]
    mass_indexing = mass_sorted[order]
    intensity_indexing = intensity_sorted[order]
    eic_index_center = np.zeros(len(mass_indexing))
    eic_index_offset = np.zeros(len(mass_indexing))
    index_counter = 0

    for idx in tqdm(range(0, len(mass_indexing[intensity_indexing>intensity_threshold]))):
        seed_mass, mass_error, left_idx, right_idx=get_seed_mass(mass_indexing[idx],mass_sorted )
        if is_in_ranges(seed_mass, eic_index_center[0:index_counter], eic_index_offset[0:index_counter])== False:

            ion_trace = flash_eic(seed_mass, mass_sorted, intensity_sorted, index_sorted,
                                  mass_error=mass_error)
            eic_index_center[index_counter]=seed_mass
            eic_index_offset[index_counter]=mass_error
            index_counter=index_counter+1
            try:
                peaks_all, raw_apex_idx_all,reci_snrs_all  = detect_all_peaks(ion_trace, n_neighbor=n_neighbor,intensity_threshold=intensity_threshold)
            except:
                print(seed_mass, mass_error)
                return
            if len(peaks_all)>0:
                for p, r, a in zip(peaks_all, reci_snrs_all, raw_apex_idx_all):
                    pmz_statistics = guess_pmz(seed_mass, mass_sorted,
                                               intensity_sorted, index_sorted, left_idx, right_idx, int(a), mass_error)

                    if pmz_statistics[0]==pmz_statistics[0]:
                        pmz[counter]=pmz_statistics[0]
                        anchor_pmzs[counter]=pmz_statistics[2]
                        apex_intensity[counter]=pmz_statistics[1]
                        ms1_idx[counter]=a
                        reci_snrs[counter]=r
                        rt_start[counter]=rt_list[p[0]]
                        rt_end[counter]=rt_list[p[2]]
                        n_scans[counter]=p[2]-p[0]-1
                        peak_range[counter]=[p[0],a,p[2]]
                        eic_center[counter] = eic_center[counter] + seed_mass
                        eic_offset[counter]=eic_offset[counter]+mass_error
                        try:
                            rt[counter]=gaussian_estimator(tuple([int(a)-1,int(a), int(a)+1]),rt_list, ion_trace)
                            if rt[counter]!=rt[counter]:
                                rt[counter]=rt_list[int(a)]
                        except:
                            rt[counter]=rt_list[int(a)]

                        counter = counter + 1

    eic_center = eic_center[0:counter]
    eic_offset = eic_offset[0:counter]
    pmz = pmz[0:counter]
    rt = rt[0:counter]
    rt_start = rt_start[0:counter]
    rt_end = rt_end[0:counter]
    apex_intensity = apex_intensity[0:counter]
    n_scans = n_scans[0:counter]
    peak_range = peak_range[0:counter]
    reci_snrs=reci_snrs[0:counter]
    anchor_pmzs=anchor_pmzs[0:counter]
    df = pd.DataFrame(zip(pmz,anchor_pmzs, eic_center, eic_offset,
                          rt, rt_start, rt_end,
                          apex_intensity,
                          n_scans,peak_range,reci_snrs),
                      columns=['precursor_mz','anchor_pmz','eic_center', 'eic_offset',
                               'rt_apex', 'rt_start', 'rt_end',
                               'ms1_intensity',
                               'n_scnas', 'ms1_scan_range', 'reci_snr'])
    return(df)

def get_current_peak(peaks_all, all_apex_intensity,intensity_list_smoothed):
    current_peak_idx = np.argmax(all_apex_intensity)
    apex_peak = np.array(peaks_all[current_peak_idx])
    apex_peak_range = [current_peak_idx,current_peak_idx]
    l = 1
    r = 1
    while current_peak_idx-l>0:
        left_peak = peaks_all[current_peak_idx-l]
        if apex_peak[0]-left_peak[2]<=2 and intensity_list_smoothed[apex_peak[0]]>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[left_peak[1]])/1.3:
            apex_peak[0]=left_peak[0]
            apex_peak_range[0]=current_peak_idx-l
            l = l+1
        else:
            break
    while current_peak_idx+r<len(peaks_all):
        right_peak = peaks_all[current_peak_idx+r]
        if right_peak[0]-apex_peak[2]<=2 and intensity_list_smoothed[apex_peak[2]]>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[right_peak[1]])/1.3:
            apex_peak[2]=right_peak[2]
            apex_peak_range[1]=current_peak_idx+r
            r = r+1
        else:
            break
    peaks_all = peaks_all[0:apex_peak_range[0]]+peaks_all[apex_peak_range[1]+1:]
    all_apex_intensity = np.array([intensity_list_smoothed[x[1]]for x in peaks_all])
    return apex_peak, peaks_all, all_apex_intensity
def detect_all_peaks(intensity_list, intensity_threshold = 30000, n_neighbor=2, return_all = False, return_baseline = False):
    intensity_list_smoothed = np.array(moving_average(intensity_list, n= n_neighbor))
    peaks_all =get_peaks(intensity_list_smoothed)
    if len(peaks_all)==0:
        return([], np.NAN, np.NAN)
    peaks_return = [None]*len(peaks_all)
    reci_snrs = np.zeros(len(peaks_all))
    raw_apex_idx = np.zeros(len(peaks_all), dtype=int)
    # smoothed_apex_intensity =np.zeros(len(peaks_all))
    idx = 0
    all_apex_intensity = np.array([intensity_list_smoothed[x[1]]for x in peaks_all])
    if np.max(all_apex_intensity)<=intensity_threshold:
        return ([], np.NAN, np.NAN)
    while np.max(all_apex_intensity)>intensity_threshold:
        apex_peak, peaks_all, all_apex_intensity= get_current_peak(peaks_all,all_apex_intensity, intensity_list_smoothed )
        if idx == 0:
            fwhm = find_peak_width(apex_peak, intensity_list_smoothed, ratio = 1/1.5)
            fwhm = np.ceil(fwhm*1.32)
            # print(fwhm)
        peaks_return[idx]=(apex_peak)
        # smoothed_apex_intensity[idx]=intensity_list_smoothed[apex_peak[1]]
        raw_apex_idx[idx]=(apex_peak[0]+np.argmax(intensity_list[apex_peak[0]:apex_peak[2]+1]))
        idx = idx+1
        if len(peaks_all)==0:
            break
    baseline = pybaselines.smooth.snip(intensity_list,  decreasing = True, max_half_windowint = fwhm, smooth_half_window = int(fwhm))[0]
    # print(fwhm)
    # baseline = pybaselines.classification.dietrich(intensity_list_smoothed, smooth_half_window =0,poly_order=2,interp_half_window=fwhm)[0]
    # baseline = pybaselines.classification.fastchrom(intensity_list_smoothed, half_window = fwhm)[0]

    baseline[baseline<0]=0
    for i in range(0, idx):
        reci_snrs[i]=np.median(baseline[peaks_return[i][0]:peaks_return[i][2]+1])/intensity_list_smoothed[peaks_return[i][1]]
    raw_apex_idx = raw_apex_idx[0:idx]

    peaks_return = peaks_return[0:idx]
    reci_snrs=reci_snrs[0:idx]
    if return_all==True:
        toggle = raw_apex_idx>-1
    else:
        toggle = reci_snrs<1/5
    if return_baseline==False:
        return np.array(peaks_return)[toggle], raw_apex_idx[toggle], reci_snrs[toggle]
    else:
        return np.array(peaks_return)[toggle], raw_apex_idx[toggle], reci_snrs[toggle], baseline
    # n_scans=n_scans[0:idx]

# def detect_all_peaks(intensity_list, intensity_threshold = 30000, n_neighbor=2, return_all = False, return_baseline = False):
#     intensity_list_smoothed = np.array(moving_average(intensity_list, n= n_neighbor))
#     peaks_all =get_peaks(intensity_list_smoothed)
#     if len(peaks_all)==0:
#         return ([],np.NAN, np.NAN)
#     peaks_return = [None]*len(peaks_all)
#     reci_snrs = np.zeros(len(peaks_all))
#     # smoothed_apex_intensity =np.zeros(len(peaks_all))
#     idx = 0
#     raw_apex_idx = np.zeros(len(peaks_all), dtype=int)
#     all_apex_intensity = np.array([intensity_list_smoothed[x[1]]for x in peaks_all])
#
#     while np.max(all_apex_intensity)>intensity_threshold:
#         current_peak_idx = np.argmax(all_apex_intensity)
#         apex_peak = np.array(peaks_all[current_peak_idx])
#         apex_peak_range = [current_peak_idx,current_peak_idx]
#         l = 1
#         r = 1
#         while current_peak_idx-l>0:
#             left_peak = peaks_all[current_peak_idx-l]
#             if apex_peak[0]-left_peak[2]<=2 and intensity_list_smoothed[apex_peak[0]]>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[left_peak[1]])/1.3:
#                 apex_peak[0]=left_peak[0]
#                 apex_peak_range[0]=current_peak_idx-l
#                 l = l+1
#             else:
#                 break
#         while current_peak_idx+r<len(peaks_all):
#             right_peak = peaks_all[current_peak_idx+r]
#             if right_peak[0]-apex_peak[2]<=2 and intensity_list_smoothed[apex_peak[2]]>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[right_peak[1]])/1.3:
#                 apex_peak[2]=right_peak[2]
#                 apex_peak_range[1]=current_peak_idx+r
#                 r = r+1
#             else:
#                 break
#         peaks_return[idx]=(apex_peak)
#         # smoothed_apex_intensity[idx]=intensity_list_smoothed[apex_peak[1]]
#         raw_apex_idx[idx]=(apex_peak[0]+np.argmax(intensity_list[apex_peak[0]:apex_peak[2]+1]))
#         if idx ==0:
#             fwhm = find_single_fwhm(apex_peak, intensity_list_smoothed)
#
#         idx = idx+1
#         # all_apex_intensity = np.concatenate((all_apex_intensity[0:apex_peak_range[0]],all_apex_intensity[apex_peak_range[1]+1:]), axis=None)
#         peaks_all = peaks_all[0:apex_peak_range[0]]+peaks_all[apex_peak_range[1]+1:]
#         all_apex_intensity = np.array([intensity_list_smoothed[x[1]]for x in peaks_all])
#
#         if len(peaks_all)==0:
#             break
#     if idx == 0:
#         return([], np.NAN, np.NAN)
#     # n_scans=n_scans[0:idx]
#     baseline = pybaselines.smooth.snip(intensity_list_smoothed, max_half_window = fwhm, decreasing = True)[0]
#
#     # baseline = pybaselines.classification.dietrich(intensity_list_smoothed)[0]
#     baseline = np.array([0 if x <0 else x for x in baseline])
#     for i in range(0, idx):
#         reci_snrs[i]=np.median(baseline[peaks_return[i][0]:peaks_return[i][2]+1])/intensity_list_smoothed[peaks_return[i][1]]
#     peaks_return=peaks_return[0:idx]
#     reci_snrs = reci_snrs[0:idx]
#     raw_apex_idx=raw_apex_idx[0:idx]
#
#     if return_all==True:
#         return (peaks_return, reci_snrs, raw_apex_idx)
#     high_quality_idx = reci_snrs<0.2
#     peaks_return=np.array(peaks_return)
#
#     peaks_return = peaks_return[high_quality_idx]
#     reci_snrs=reci_snrs[high_quality_idx]
#     raw_apex_idx = raw_apex_idx[high_quality_idx]
#
#     if return_baseline ==True:
#         return(peaks_return,raw_apex_idx,reci_snrs, baseline )
#     else:
#         return(peaks_return,raw_apex_idx,reci_snrs )
# def detect_all_peaks(intensity_list, intensity_threshold = 30000, n_neighbor=2, split_level = 1.5, return_all = False):
#     intensity_list_smoothed = np.array(moving_average(intensity_list, n= n_neighbor))
#     # intensity_list_smoothed = intensity_list_smoothed+1
#     peaks_all =get_peaks(intensity_list_smoothed)
#     if len(peaks_all)==0:
#         return(np.NAN)
#     peaks_return = []
#     raw_apex_idx = []
#     reci_snr = []
#     baseline = get_baseline(intensity_list_smoothed)
#
#     while peaks_all == peaks_all and len(peaks_all)>0:
#         apex_peak, peaks_all = connect_one_peak(peaks_all, intensity_list_smoothed, intensity_threshold = intensity_threshold, split_level = split_level)
#         if apex_peak == apex_peak:
#             peaks_return.append(apex_peak)
#             raw_apex_idx_temp = apex_peak[0]+np.argmax(intensity_list[apex_peak[0]:apex_peak[2]+1])
#             raw_apex_idx.append(raw_apex_idx_temp)
#             reci_snr.append(np.median(baseline[apex_peak[0]:(apex_peak[2]+1)])/intensity_list[raw_apex_idx_temp])
#     peaks_return=np.array(peaks_return)
#     raw_apex_idx = np.array(raw_apex_idx)
#     reci_snr = np.array(reci_snr)
#     if return_all == True:
#         order = raw_apex_idx>-1
#     else:
#         order = reci_snr<1/3
#     return(peaks_return[order] , raw_apex_idx[order], reci_snr[order])
    # return(peaks_return[order] , raw_apex_idx[order], reci_snr[order], baseline)
# def get_baseline(intensity_list_smoothed):
#     # index_start, index_end = trim_zeros_and_get_indices(intensity_list_smoothed)
#     # has_number = intensity_list_smoothed[index_start:index_end+1]
#     # leading_zeros = intensity_list_smoothed[0:index_start]
#     # tailing_zeros = intensity_list_smoothed[index_end+1:]
#     # baseline = pybaselines.classification.dietrich(intensity_list_smoothed)[0]
#     #
#     # baseline_center[baseline_center<0]=0
#     baseline_center = np.repeat(np.median(intensity_list_smoothed), len(intensity_list_smoothed))
#     # baseline = np.concatenate([leading_zeros, baseline_center, tailing_zeros])
#     return baseline_center
def trim_zeros_and_get_indices(arr):
    """
    Remove all consecutive zeros from the start and end of a NumPy array and return the start and end indices of the remaining array.

    :param arr: Input NumPy array.
    :return: A tuple (start_index, end_index) of the remaining array.
    """
    # Trim zeros from both ends
    trimmed_arr = np.trim_zeros(arr, 'fb')

    # Find the start and end indices of the non-zero part in the original array
    if trimmed_arr.size > 0:
        start_index = np.where(arr == trimmed_arr[0])[0][0]
        end_index = np.where(arr == trimmed_arr[-1])[0][-1]
        return start_index, end_index
    else:
        # If the trimmed array is empty, return None as indices
        return None, None
def find_peak_width(peak, intensity_list_smoothed, ratio = 0.5):
    peak = peak[1]
    peak_height = intensity_list_smoothed[peak]
    half_max = peak_height *ratio
    left_idx = np.where(intensity_list_smoothed[:peak] <= half_max)[0]
    if len(left_idx) > 0:
        left_idx = left_idx[-1]
    else:
        left_idx = peak
    right_idx = np.where(intensity_list_smoothed[peak:] <= half_max)[0]
    if len(right_idx) > 0:
        right_idx = peak + right_idx[0]
    else:
        right_idx = peak
    fwhm_width = right_idx - left_idx
    return(fwhm_width)
# def connect_one_peak(peaks_all, intensity_list_smoothed,split_level = 1.5, intensity_threshold = 30000):
#     all_apex_idx = np.array([arr[1] for arr in peaks_all])
#     all_apex_intensity = intensity_list_smoothed[all_apex_idx]
#     if np.max(all_apex_intensity)<intensity_threshold:
#         return (np.NAN, np.NAN)
#     current_peak_idx = np.argmax(all_apex_intensity)
#     l = 1
#     r = 1
#     # searching left
#     apex_peak = (peaks_all[current_peak_idx])
#     apex_peak_range = [current_peak_idx,current_peak_idx]
#
#     while current_peak_idx-l>=0:
#         left_peak = peaks_all[current_peak_idx-l]
#         if apex_peak[0]-left_peak[2]<=2 and max(intensity_list_smoothed[apex_peak[0]], intensity_list_smoothed[left_peak[2]])>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[left_peak[1]])/split_level:
#             apex_peak[0]=left_peak[0]
#             apex_peak_range[0]=current_peak_idx-l
#             l = l+1
#         else:
#             break
#     while current_peak_idx+r<=len(peaks_all)-1:
#         right_peak = peaks_all[current_peak_idx+r]
#         if right_peak[0]-apex_peak[2]<=2 and max(intensity_list_smoothed[apex_peak[2]], intensity_list_smoothed[right_peak[0]])>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[right_peak[1]])/split_level:
#             apex_peak[2]=right_peak[2]
#             apex_peak_range[1]=current_peak_idx+r
#             r = r+1
#         else:
#             break
#     peaks_all_return = peaks_all[0:apex_peak_range[0]]+peaks_all[apex_peak_range[1]+1:]
#     return(apex_peak, peaks_all_return)






def guess_pmz(target_mass,mass_sorted, intensity_sorted,index_sorted,idx_start, idx_end, peak_apex, mass_error,guess_step = 1):
    # idx_start, idx_end = mass_sorted.searchsorted([seed_mass-mass_tolerance,seed_mass+mass_tolerance ])
    mass_range = mass_sorted[idx_start:idx_end]
    index_range = index_sorted[idx_start:idx_end]
    intensity_range = intensity_sorted[idx_start:idx_end]
    if np.max(intensity_range)==0:
        return(np.NAN, np.NAN)
    pmz_candidates = np.zeros(2*guess_step+1)
    intensity_candidates = np.zeros(2*guess_step+1)
    pmz_idx = 0
    for i in range(peak_apex-guess_step, peak_apex+guess_step+1):
        if len(intensity_range[index_range==i])>0:
            difference_array = np.absolute(mass_range[index_range==i]-target_mass)

            mass_anchor = difference_array.argmin()
            pmz_candidates[pmz_idx]=mass_range[index_range==i][mass_anchor]
            if i == peak_apex:
                anchor_pmz = mass_range[index_range==i][mass_anchor]
            intensity_candidates[pmz_idx]=intensity_range[index_range==i][mass_anchor]
            pmz_idx = pmz_idx+1
        # if i == peak_apex:
        #     apex_intensity = intensity_candidates[pmz_idx]

    if pmz_idx<(2*guess_step+1):
        return(np.NAN, np.NAN,np.NAN)
    pmz_candidates=pmz_candidates[0:pmz_idx]
    intensity_candidates=intensity_candidates[0:pmz_idx]
    weighted_pmz = np.sum([x*y/np.sum(intensity_candidates) for x, y in zip(pmz_candidates, intensity_candidates)])
    apex_intensity = flash_eic(weighted_pmz, mass_sorted, intensity_sorted, index_sorted, mass_error = mass_error)[peak_apex]
    return (weighted_pmz, apex_intensity, anchor_pmz)

def map_ms2(features_df, ms2):
    msms = []
    msms_most_fragmented = []
    ms2_scan_idx = []
    rt_offset = []
    rt_offset_most_fragmented = []
    ms2_pmz = []
    ms2_pmz_intensity = []
    ms2_pmz_most_fragmented = []
    ms2_pmz_intensity_most_fragmented = []
    ms2_working = ms2.copy()
    features_df.sort_values(by = 'ms1_intensity', ascending=False, inplace=True)
    ms2_working.sort_values(by = 'precursor_mz', ascending=True, inplace=True)
    for index, row in features_df.iterrows():
        if len(ms2_working)==0:
            features_df['msms']=np.NAN
            features_df['msms_pmz']=np.NAN
            features_df['msms_pmz_intensity']=np.NAN
            features_df['msms_idx']=np.NAN
            features_df['rt_offset']=np.NAN
            features_df['msms_mf']=np.NAN

            features_df['msms_mf_pmz']=np.NAN
            features_df['msms_mf_pmz_intensity']=np.NAN
            features_df['rt_offset_mf']=np.NAN
            return features_df
            break
        pmz_matched = quick_search_sorted(ms2_working, 'precursor_mz', row['precursor_mz']-0.01, row['precursor_mz']+0.01)
        ms2_under_peak = quick_search_values(pmz_matched, 'ms1_rt', row['rt_start'], row['rt_end'])

        if len(ms2_under_peak)>0:

            ms2_under_peak['intensity_coef']=ms2_under_peak['ms1_precursor_intensity']*ms2_under_peak['peak_purity']
            ms2_under_peak.sort_values(by = 'intensity_coef', ascending=False, inplace = True)
            fragment_intensity = [so.get_fragment_intensity(row['peaks'],row['precursor_mz'] ) for index, row in ms2_under_peak.iterrows()]
            most_fragmented_row = ms2_under_peak.iloc[np.argmax(fragment_intensity)]
            msms_most_fragmented.append(most_fragmented_row['peaks'])
            msms.append((ms2_under_peak.iloc[0]['peaks']))
            ms2_pmz_intensity.append(ms2_under_peak.iloc[0]['precursor_intensity'])
            rt_offset.append(abs(row['rt_apex']-ms2_under_peak.iloc[0]['rt']))
            ms2_pmz.append(ms2_under_peak.iloc[0]['precursor_mz'])
            ms2_pmz_most_fragmented.append(most_fragmented_row['precursor_mz'])
            rt_offset_most_fragmented.append(abs(row['rt_apex']-most_fragmented_row['rt']))
            ms2_pmz_intensity_most_fragmented.append(most_fragmented_row['precursor_intensity'])
            ms2_scan_idx.append(int(ms2_under_peak.iloc[0]['scan_idx']))
            ms2_working.drop(ms2_under_peak.index.tolist(), inplace=True)
        else:
            msms.append(np.NAN)
            ms2_pmz.append(np.NAN)
            ms2_pmz_intensity.append(np.NAN)
            ms2_pmz_intensity_most_fragmented.append(np.NAN)
            ms2_pmz_most_fragmented.append(np.NAN)
            ms2_scan_idx.append(np.NAN)
            rt_offset.append(np.NAN)
            msms_most_fragmented.append(np.NAN)
            rt_offset_most_fragmented.append(np.NAN)
    features_df['msms']=msms
    features_df['msms_pmz']=ms2_pmz
    features_df['msms_pmz_intensity']=ms2_pmz_intensity
    features_df['msms_idx']=ms2_scan_idx
    features_df['rt_offset']=rt_offset
    features_df['msms_mf']=msms_most_fragmented

    features_df['msms_mf_pmz']=ms2_pmz_most_fragmented
    features_df['msms_mf_pmz_intensity']=ms2_pmz_intensity_most_fragmented
    features_df['rt_offset_mf']=rt_offset_most_fragmented

    return(features_df)


import itertools
def build_index(ms1, use_binned = False):
    # ms1.reset_index(inplace = True, drop = True)
    if use_binned == False:
        col = 'peaks'
    else:
        col = 'peaks_binned'
    mass_nested = [None]*len(ms1)
    intensity_nested = [None]*len(ms1)
    index_nested = [None]*len(ms1)
    rt_list = np.zeros(len(ms1))
    for index, row in (ms1.iterrows()):
        mass_temp, intensity_temp = row[col].T
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
    # loc_sorted = np.arange(len(mass_sorted))
    return(mass_sorted, intensity_sorted, index_sorted, rt_list)
def flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted, mass_error=0.005, gap_fill = True):
    index_start, index_end = mass_sorted.searchsorted([pmz-mass_error, pmz+mass_error+1E-9])
    index_range = index_sorted[index_start:index_end]

    intensity_range = intensity_sorted[index_start:index_end]

    intensity_list = np.zeros(np.max(index_sorted)+1)
    for idx in range(0,len(index_range)):
        intensity_list[index_range[idx]]= intensity_list[index_range[idx]]+intensity_range[idx]
    if gap_fill == True:
        intensity_list = gap_filling(intensity_list, max_gap=2)
    return(intensity_list)

# def flash_eic_single(pmz, mass_sorted, intensity_sorted, index_sorted, mass_error=0.005):
#     index_start, index_end = mass_sorted.searchsorted([pmz-mass_error, pmz+mass_error+1E-9])
#     intensity_range = intensity_sorted[index_start:index_end]
#     index_range = index_sorted[index_start:index_end]
#     mass_range = mass_sorted[index_start:index_end]
#
#     mass_offset = np.abs(mass_range-pmz)
#     order = np.argsort(mass_offset)
#     index_range = index_range[order]
#     intensity_range = intensity_range[order]
#     intensity_list = np.zeros(np.max(index_sorted)+1)
#     counter = 0
#     for idx in range(0,len(index_range)):
#         if intensity_list[index_range[idx]]==0:
#             intensity_list[index_range[idx]]= intensity_list[index_range[idx]]+intensity_range[idx]
#             counter = counter+1
#     intensity_list = gap_filling(intensity_list)
#     return(intensity_list)
def moving_average( x, n=2):
    window = np.ones(2*n + 1) / (2*n + 1)

    # Use 'same' mode in convolve to ensure the output array has the same size as the input array
    # 'valid' mode would reduce the size of the array to avoid edge effects
    # 'full' mode would increase the size to include all possible overlaps
    smoothed_arr = np.convolve(x, window, mode='same')

    return smoothed_arr
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
    # gradient = np.diff(intensity_list)
    left_edge_idx = cur_apex_idx-1
    right_edge_idx = cur_apex_idx+1
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

    return([left_edge_idx, cur_apex_idx, right_edge_idx])


#%%
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
def gap_filling(intensity_list, max_gap = 2):
    zero_ranges = zero_runs(intensity_list)
    intensity_list_twisted = intensity_list.copy()
    for zr in zero_ranges:
        if zr[1]-zr[0]<=max_gap and zr[0]!=0 and zr[1]!=len(intensity_list):
            gradient = (intensity_list[zr[1]]-intensity_list[zr[0]-1])/(zr[1]-zr[0]+1)
            for j in range(0, zr[1]-zr[0]):

                intensity_list_twisted[j+zr[0]]=intensity_list_twisted[j+zr[0]]+ intensity_list[zr[0]-1]+gradient*(j+1)
    return(intensity_list_twisted)
def zero_runs(a):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges
# def find_roi(ion_trace):
#     apex_index = np.argmax(ion_trace)
#     ion_trace_smoothed = trx.moving_average(ion_trace, 10)
#     if np.min(ion_trace_smoothed)>0:
#         return([0, apex_index, len(ion_trace)-1])
#     apex_smoothed = np.argmax(ion_trace)
#     left_side = ion_trace_smoothed[0:apex_smoothed]
#     right_side = ion_trace_smoothed[apex_smoothed+1:]
#     if np.min(left_side)>0:
#         left_bound = 0
#     else:
#         # for i in range(1, apex_smoothed+1):
#         #     if ion_trace_smoothed[apex_smoothed-i]==0:
#         #         left_bound= apex_smoothed-i
#         #         break
#         zero_position_left = np.where(left_side == 0)[0]
#         # index = (apex_smoothed-zero_position_left).argmin()
#         left_bound = zero_position_left[-1]
#     if np.max(right_side)>0:
#         right_bound = len(ion_trace)-1
#     else:
#         # for l in range(1, len(apex_smoothed)-apex_index):
#         #     if ion_trace_smoothed[apex_smoothed+l]==0:
#         #         right_bound = apex_smoothed+l
#         #         break
#         zero_position_right = np.where(right_side == 0)[0]
#         # index = (zero_position_right-apex_smoothed).argmin()
#         right_bound = zero_position_right[0]
#     return(np.array([left_bound, apex_index, right_bound]))
