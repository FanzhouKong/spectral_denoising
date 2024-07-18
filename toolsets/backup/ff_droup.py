import time

if __name__ == '__main__':
    freeze_support()
import numpy as np
# import pymzml
import pandas as pd
import gc


import os
import toolsets.spectra_operations as so
from toolsets.search import string_search, quick_search_values, num_search
import toolsets.helpers as helpers
import toolsets.file_io as io
from itertools import repeat
import pymzml
from scipy.signal import find_peaks
import toolsets.parallel_functions as pf
# from toolsets.parallel_functions import _extract_mzml_info
from multiprocessing import Pool, freeze_support
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
import contextlib
import concurrent.futures
from matplotlib import rcParams


from concurrent import futures
import matplotlib.pyplot as plt

import seaborn as sns

import plotly.express as px

def feature_finding(ms1, ms2):
    bins = get_mz_bin(ms2)
    with Pool() as pool:
        results= pool.starmap(get_feature, zip(repeat(ms2), bins, repeat(ms1)))
    # print('i exited pool okay')
    features = pd.DataFrame()
    for result in results:
        features = pd.concat([features, result], ignore_index=True)
    # features.drop_duplicates(subset=['ms2_range_idx'], inplace=True)
    features.reset_index(inplace=True, drop=True)
    return(features)
    # for feature

def get_feature(ms2,current_pmz_bin,ms1):
    current_cluster = quick_search_values(ms2, 'precursor_mz', value_start=current_pmz_bin-0.005, value_end=current_pmz_bin+0.005)
    features_all = pd.DataFrame()
    rt_list, intensity_list = get_EIC_list(ms1, current_pmz_bin)
    while len(current_cluster)>0:
       # print('length of file is: ', len(current_cluster))
       seed_ms2 = current_cluster.iloc[np.argmax(current_cluster['ms1_precursor_intensity'])]
       # rt_list, intensity_list = get_EIC_list(ms1, seed_ms2['ms1_pmz'])
       # rt_cnt, int_list = ff.centroid_data(rt_list, intensity_list)

       peak_list = get_peaks(intensity_list)
       target_peak_idx = find_most_close(peak_list, rt_list, seed_ms2['rt'], return_index=True)


       if len(peak_list) >0 and target_peak_idx != -1:
           target_peak = connect_peaks(peak_list, target_peak_idx, intensity_list, rt_list)
           rt_start, rt_end = rt_list[target_peak[0]], rt_list[target_peak[2]]
           rt_apex, int_apex = get_centroid(target_peak, rt_list, intensity_list)
           current_cluster_rt_sorted = current_cluster.sort_values(by = 'rt', ascending=True)
           ms2_under_peaks  = quick_search_values(current_cluster_rt_sorted, 'rt', value_start=rt_start, value_end=rt_end)

           ms2_under_peaks['rt']=rt_apex
           ms2_under_peaks['rt_offset']=abs(rt_apex-seed_ms2['rt'])
           ms2_under_peaks['rt_start']=rt_start
           ms2_under_peaks['rt_end']=rt_end
           ms2_under_peaks['precursor_mz']=(ms2_under_peaks['precursor_mz']*ms2_under_peaks['ms1_precursor_intensity']).sum()/(ms2_under_peaks['ms1_precursor_intensity'].sum())
           # pmzs = []
           # intensities = []
           # for i in range(target_peak[0], target_peak[2]+1):
           #     pmz_ms1, precursor_intensity, mz_offset= _extract_precursor_intensity(ms1.iloc[i]['peaks'], current_pmz_bin)
           #     pmzs.append(pmz_ms1)
           #     intensities.append(precursor_intensity)
           # ms2_under_peaks['precursor_mz']=sum(x * y for x, y in zip(pmzs, intensities))/sum(intensities)


           ms2_under_peaks['peak_apex_intensity']=int_apex

           ms2_under_peaks['peaks']=so.weighted_average_spectra(ms2_under_peaks, weight_col= 'ms1_precursor_intensity')
           ms2_under_peaks['peak_range_idx'] = ms2_under_peaks.apply(lambda x: list(target_peak), axis=1)
           ms2_under_peaks['ms2_range_idx'] = ms2_under_peaks.apply(lambda x: list(ms2_under_peaks['scan_idx'].unique()), axis=1)
           ms2_under_peaks['pmz_bin']=current_pmz_bin
           ms2_under_peaks.drop(['scan_idx', 'cycle', 'isolation_window', 'ms1_pmz', 'ms1_precursor_intensity', 'peak_purity', 'mz_offset', 'ms1_pmz', 'ms1_rt'], axis =1, inplace = True)
           current_cluster.drop(ms2_under_peaks.index, inplace = True)
           features_all=pd.concat([features_all,pd.DataFrame([ms2_under_peaks.iloc[0]]) ], axis=0)
       else:
           # print('npmax', np.argmax(current_cluster['ms1_precursor_intensity']))
           current_cluster = string_search(current_cluster, 'scan_idx', seed_ms2.scan_idx, reverse=True)
           # break
    return(features_all)
    # if len(features_all)>0:

    #     features_tidy = pd.DataFrame()

    #     for rt in features_all['rt'].unique():
    #        df_temp = string_search(features_all, 'rt', rt)
    #        df_temp.sort_values(by = 'rt_offset', inplace = True, ascending=True)
    #        features_tidy = features_tidy.append(df_temp.iloc[0], ignore_index=True)
    #     return(features_tidy)
    # else:
    #     return(features_all)

   

   
       




   


def connect_peaks(peak_list, target_peak_idx, intensity_list, rt_list):
    apex_peak = list(peak_list[target_peak_idx])
    # peak_start = apex_peak[0]
    # peak_apex = apex_peak[1]
    # peak_end = apex_peak[2]

    for i in range(1, target_peak_idx+1):
       left_peak = peak_list[target_peak_idx-i]
       # break
       if apex_peak[0]-left_peak[2]<=2 and intensity_list[apex_peak[0]]>min(intensity_list[apex_peak[1]], intensity_list[left_peak[1]])/1.3 and rt_list[apex_peak[2]]-rt_list[apex_peak[0]]<40/60:
          apex_peak[0] = left_peak[0]
          apex_peak[1]=np.argmax(intensity_list[apex_peak[0]:apex_peak[2]])+apex_peak[0]
       else:
          break

    for i in range(1, len(peak_list)-target_peak_idx):
       right_peak = peak_list[target_peak_idx+i]
       # break
       if right_peak[0]-apex_peak[2]<=2 and intensity_list[apex_peak[2]]>min(intensity_list[apex_peak[1]], intensity_list[right_peak[1]])/1.3 and rt_list[apex_peak[2]]-rt_list[apex_peak[0]]<40/60:
          apex_peak[2] = right_peak[2]
          apex_peak[1]=np.argmax(intensity_list[apex_peak[0]:apex_peak[2]])+apex_peak[0]
       else:
          break
    # apex_peak[1]=np.argmax(intensity_list[apex_peak[0]:apex_peak[2]])+apex_peak[0]
    return(tuple(apex_peak))
def process_mzml(mzml_path, parent_dir =  None, rt_max = 20,if_mix = True, with_ms1 = True):
    if if_mix == True and parent_dir != None:
        mix = mzml_path
        mzml_base_name = helpers.find_files(parent_dir, mix)
        mzml_path = os.path.join(parent_dir, mzml_base_name)
    if if_mix == False and parent_dir!= None:
        mzml_path = os.path.join(parent_dir,mzml_path)
        # print(mzml_path)
        # print('i finished loading mzml file')
    try:
        ms1_2 = load_mzml_data(mzml_path, rt_max = rt_max)
    except:
        print('the mzml file or mix name you passed might be wrong')
        # raise Error
    ms2 = string_search(ms1_2, 'ms_level', 2)
    
    
    
    ms1 = string_search(ms1_2, 'ms_level', 1)
    ms1.sort_values(by='rt', inplace=True)
    ms1.reset_index(inplace=True, drop=True)
    ms1_intensity = []
    ms1_rt = []
    peak_purity = []
    ms1_pmz=[]
    mz_offsets = []
    for index, row in ms2.iterrows():
        ms1_scan = string_search(ms1, 'cycle', row['cycle']).iloc[0]
        pmz_ms1, precursor_intensity, mz_offset= _extract_precursor_intensity(ms1_scan['peaks'], row['precursor_mz'])
        ms1_rt.append(ms1_scan['rt'])

        # EIC_rt_temp, EIC_int_temp = get_EIC_list(ms1, pmz_ms1, step = 0.005)
        # EIC_rt.append(EIC_rt_temp)
        # EIC_int.append(EIC_int_temp)
        isolation_window_intensity = _extract_ms1_intensity(ms1_scan['peaks'], row['isolation_window'][0], row['isolation_window'][1])
        if isolation_window_intensity!= 0:
            peak_purity.append(precursor_intensity/isolation_window_intensity)
        else:
            peak_purity.append(0)
        ms1_intensity.append(precursor_intensity)

        ms1_pmz.append(pmz_ms1)
        mz_offsets.append(mz_offset)
    # print('i am done calculating ms1 intensity and peak purity')
    # ms2['EIC_rt']=EIC_rt
    # ms2['EIC_int']=EIC_int
    ms2['ms1_pmz']=ms1_pmz
    ms2['ms1_rt']=ms1_rt
    ms2['ms1_precursor_intensity'] =ms1_intensity
    ms2['peak_purity']=peak_purity
    ms2['mz_offset']=mz_offsets
    # ms2['mix']=mix
    if if_mix== True:
        ms2['base_name']=mzml_base_name
        ms2['mix']=mix
    else:
        ms2['mix']= os.path.basename(mzml_path)
        ms2['base_name']=os.path.basename(mzml_path)
    ms2.sort_values(by = ['ms1_pmz', 'ms1_precursor_intensity'], inplace = True)
    ms2.reset_index(inplace=True, drop=True)
    if with_ms1 == True:
        return(ms1, ms2)
    else:
        return(ms2 )


def get_EIC_list(ms1, pmz, step = 0.005):
    rt_list = []
    intensity_list = []

    for index, row in ms1.iterrows():
        rt_list.append(row['rt'])
        intensity_list.append(_extract_ms1_intensity(row['peaks'], mz_lower=pmz-step, mz_upper=pmz+step))
    return(rt_list, intensity_list)
def auto_EIC(mix, parent_dir,pmz, vlines_location_1=[], vlines_location_2=[] , rt_start = -1, rt_end = -1, rt_max = 5):
    ms1, ms2 = process_mzml(mix, parent_dir, if_mix=True, with_ms1=True, rt_max = rt_max)
    rt_list, intensity_list = get_EIC_list(ms1, pmz)
    EIC(rt_list, intensity_list, vlines_location_1=vlines_location_1,vlines_location_2 = vlines_location_2, rt_start=rt_start, rt_end=rt_end)
def EIC(rt_list, intensity_list,
    # parent_dir= None, if_mix = False, 
    rt_start = -1, rt_end = -1, adjusted_height = -1, vlines_location_1 = [], vlines_location_2 = [],
        savepath = None):
    fig, ax = plt.subplots(
    figsize = (9, 5)
                      )
    ax= sns.lineplot(x = rt_list, y = intensity_list)
    if rt_start != -1 and rt_end != -1:
        index_start = np.searchsorted(rt_list, rt_start,side = 'left')
        index_end = np.searchsorted(rt_list, rt_end,side = 'right')
        # rt_list = rt_list[index_start:index_end+1]

        ax.set_xlim(rt_start, rt_end)
        ax.set_ylim(0, adjusted_height)
        ax.set_ylim(0, np.max(intensity_list[index_start:index_end])*1.1)
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
    # ax.set(xticklabels=[], yticklabels = [])
    fig.tight_layout()
    # ax.patch.set_linewidth(1)
    # ax.set_edgecolor("black")
    if savepath != None:
        plt.savefig(savepath, dpi = 300,facecolor = 'white', edgecolor = 'black')
    # return(rt_list, intensity_list)
    # ms1_intensity = []
    # peak_purity = []



def load_mzml_data(file: str, n_most_abundant=400, rt_max = 5) -> tuple:
    """Load data from an mzml file as a dictionary.

    Args:
        filename (str): The name of a .mzml file.
        n_most_abundant (int): The maximum number of peaks to retain per MS2 spectrum.
        callback (callable): A function that accepts a float between 0 and 1 as progress. Defaults to None.

    Returns:
        a pandas dataframe with all the raw data

    """
    
    # np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    # mass_list = []
    # int_list = []
    ms_list = []
    prec_mzs_list = []
    mono_mzs_list = []
    charge_list = []
    select_windows_list = []
    cycles = []
    cycle = -1
    scan_idx = []
    cur_idx = 0
    peak_list = []
    rt_list = []
    id = 1
    # specs = pymzml.run.Reader(file, obo_version="4.1.33")
    # with futures.ThreadPoolExecutor() as executor:
    # # submit all tasks
    # for result in executor.map(task):
    #     print(result)
    for spec in pymzml.run.Reader(file, obo_version="4.1.33", build_index_from_scratch=True):

        rt, masses, intensities, ms_order, prec_mass, mono_mz, charge, (prec_windows_lower, prec_windows_upper)=_extract_mzml_info(spec)
        if rt_max!= None and rt>rt_max:
            break

        to_keep = intensities > 0
        masses = masses[to_keep]
        intensities = intensities[to_keep]
        if ms_order == 1:
            cycle = cycle+1
        cycles.append(cycle)
        if not np.all(masses[:-1] <= masses[1:]):
            order = np.argsort(masses)
            masses = masses[order]
            intensities = intensities[order]

        # Only keep top n_most_abundant peaks
        if ms_order == 2 and len(masses) > n_most_abundant:
            sortindex = np.argsort(intensities)[::-1][:n_most_abundant]
            sortindex.sort()
            masses, intensities = masses[sortindex], intensities[sortindex]
        peak_list.append(so.pack_spectra(masses, intensities)) 
        ms_list.append(ms_order)
        prec_mzs_list.append(prec_mass)
        mono_mzs_list.append(mono_mz)
        charge_list.append(charge)
        select_windows_list.append((prec_windows_lower, prec_windows_upper))
        scan_idx.append(cur_idx)
        rt_list.append(rt)
        cur_idx = cur_idx+1
    return(pd.DataFrame({
        'scan_idx':scan_idx,
        'cycle':cycles,
        'ms_level':ms_list,
        'precursor_mz':mono_mzs_list,
        'charge':charge_list,
        'rt':rt_list,
        'peaks':peak_list,
        'isolation_window':select_windows_list
        }))
                
    
def _calculate_mass(mono_mz: float, charge: int) -> float:
    """Calculate the precursor mass from mono mz and charge.

    Args:
        mono_mz (float): mono m/z.
        charge (int): charge.

    Returns:
        float: precursor mass.
        TBH this doesnt get a lot of use, since it assumes [M+H]+ adduct
    """
    prec_mass = mono_mz * abs(charge) - charge * 1.00727646687

    return prec_mass



def _extract_mzml_info(input_dict: dict) -> tuple:
    """Extract basic MS coordinate arrays from a dictionary.

    Args:
        input_dict (dict): A dictionary obtained by iterating over a Pyteomics mzml.read function.

    Returns:
        tuple: The rt, masses, intensities, ms_order, prec_mass, mono_mz, charge arrays retrieved from the input_dict.
            If the `ms level` in the input dict does not equal 2, the charge, mono_mz and prec_mass will be equal to 0.

    """
    rt = input_dict.scan_time_in_minutes()
    masses = input_dict.mz
    intensities = input_dict.i
    ms_order = input_dict.ms_level
    prec_mass = mono_mz = charge = 0
    if ms_order == 2 and len(input_dict.selected_precursors) > 0:
        charge = input_dict.selected_precursors[0].get("charge", 1)
        if input_dict.get_element_by_name('negative scan') is not None:
            charge = -charge
        elif input_dict.get_element_by_name('positive scan') is None:
            raise Exception("Can't determine charge")
        mono_mz = input_dict.selected_precursors[0]["mz"]
        prec_mass = _calculate_mass(mono_mz, charge)
        try:
            prec_windows_center = float(input_dict.get_element_by_name("isolation window target m/z").attrib["value"])
        except:
            prec_windows_center = mono_mz
        try:
            prec_windows_lower = prec_windows_center-float(input_dict.get_element_by_name("isolation window lower offset").attrib["value"])
            prec_windows_upper = prec_windows_center+float(input_dict.get_element_by_name("isolation window upper offset").attrib["value"])
        except:
            prec_windows_lower = prec_windows_center-0.5
            prec_windows_upper = prec_windows_center+0.5
    else:
        prec_windows_lower, prec_windows_upper = 0., 0.

    return rt, masses, intensities, ms_order, prec_mass, mono_mz, charge, (prec_windows_lower, prec_windows_upper)
def _extract_precursor_intensity(peaks, pmz):
    mass_temp, intensity_temp = so.break_spectra(peaks)
    if len(mass_temp)==0:
        return(pmz, 0, np.NAN)
    index_start = np.searchsorted(mass_temp, pmz-0.5,side = 'left')
    index_end = np.searchsorted(mass_temp, pmz+0.5,side = 'right')
    mass_temp = mass_temp[index_start:index_end]
    intensity_temp = intensity_temp[index_start:index_end]
    if len(mass_temp)==0:
        return(pmz, 0, np.NAN)
    ms1_pmz_idx = np.argmax(intensity_temp)
    ms1_pmz = mass_temp[ms1_pmz_idx]
    intensity_pmz = intensity_temp[ms1_pmz_idx]
    mz_offset = abs(ms1_pmz-pmz)
    # offsets = [abs(x-pmz) for x in mass_temp]
    # ms1_pmz_idx = np.argmin(offsets)
    # ms1_pmz = mass_temp[ms1_pmz_idx]
    # precursor_index_start = np.searchsorted(mass_temp, pmz-0.005,side = 'left')
    # precursor_index_end = np.searchsorted(mass_temp, pmz+0.005,side = 'right')
    # ms1_pmz_intensity = intensity_temp[ms1_pmz_idx]
    # mz_offset = ms1_pmz-pmz
    return(ms1_pmz, intensity_pmz, mz_offset)

def _extract_ms1_intensity(peaks, mz_lower, mz_upper):
    mass_temp, intensity_temp = so.break_spectra(peaks)
    index_start = np.searchsorted(mass_temp, mz_lower,side = 'left')
    index_end = np.searchsorted(mass_temp, mz_upper,side = 'right')
    return(np.sum(intensity_temp[index_start: index_end]))

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



import numpy as np
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
def get_edges(intensity_list, cur_apex_idx):
    left_edge_int = intensity_list[cur_apex_idx]
    # check if left of the apex is less than current apex
    for i in range(1, cur_apex_idx+1):
       if intensity_list[cur_apex_idx-i]<=left_edge_int and left_edge_int!=0:
          left_edge_int = intensity_list[cur_apex_idx-i]
          left_edge_idx = cur_apex_idx-i
       else:
          left_edge_idx = cur_apex_idx-i+1
          break
    # going right
    right_edge_int = intensity_list[cur_apex_idx]
    for i in range(1, len(intensity_list)-cur_apex_idx):
       if intensity_list[cur_apex_idx+i]<=right_edge_int and right_edge_int != 0:
          right_edge_int = intensity_list[cur_apex_idx+i]
          right_edge_idx = cur_apex_idx+i
       else:
          right_edge_idx = cur_apex_idx+i-1
          break
    return((left_edge_idx, cur_apex_idx, right_edge_idx))
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
def get_mz_bin(ms2):
    ms2_working = ms2.copy()
    bins = []
    ms2_working.sort_values(by = 'precursor_mz', inplace = True, ignore_index = True)
    while len(ms2_working)>0:
        current_cluster = quick_search_values(ms2_working, 'precursor_mz', 
            ms2_working.iloc[np.argmax(ms2_working['ms1_precursor_intensity'])]['precursor_mz']-0.005, 
            ms2_working.iloc[np.argmax(ms2_working['ms1_precursor_intensity'])]['precursor_mz']+0.005)
        ms2_working.drop(current_cluster.index, inplace=True)
        bins.append(current_cluster.iloc[np.argmax(current_cluster['ms1_precursor_intensity'])]['precursor_mz'])
    return(bins)
