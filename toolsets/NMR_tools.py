import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# from toolsets.ff_droup import get_peaks
from toolsets.search import string_search, quick_search_values
import scipy.stats
import statsmodels.api as sm
from scipy import interpolate
import toolsets.spectra_operations as so
import toolsets.spectra_plotter as sp
step = 0.000275
def nmr_processing(std_smoothed, fraction_smoothed):
    model_cluster, cluster_shift_statistics = get_model_cluster(std_smoothed, fraction_smoothed)
    fraction_peaks = get_all_peaks(fraction_smoothed)
    matched_statistics = get_matched_statistics(model_cluster, fraction_smoothed)
    ratio = matched_statistics['ratio'].min()
    # print(ratio)
    std_smoothed_c = std_smoothed.copy()
    std_smoothed_c['intensity_adjusted']=std_smoothed_c['intensity_adjusted']*ratio
    # ax = nmr_stack(fraction_smoothed, std_smoothed_c)
    cluster_shift_statistics['simulated_area']=cluster_shift_statistics['area']*ratio
    cluster_shift_statistics['simulated_intensity']=cluster_shift_statistics['intensity']*ratio
    # plt.close()
    matched = []
    for index, row in cluster_shift_statistics.iterrows():
        region_peaks = quick_search_values(fraction_peaks, 'peak_apex', row['peak_apex']+row['shift']-0.02, row['peak_apex']+row['shift']+0.02, ifsorted=False)
        # if region_peaks['area'].max()*1.2>row['simulated_area']:
        if region_peaks['intensity'].max()*1.4>row['simulated_intensity']:
            matched.append(True)
        else:
            matched.append(False)
    cluster_shift_statistics['matched']=matched
    return cluster_shift_statistics, std_smoothed_c, matched_statistics

def calculate_explained_area(results):
    return(results[results['matched']==True]['area'].sum()/results['area'].sum())
def readin_nmr_profile(path, dir=None, if_smooth = True):
    col_names = ['cs', 'intensity', 'discard']
    if dir != None:
        data = pd.read_csv(os.path.join(dir, path), sep= '\t', header=None,names=col_names, index_col=None)
    else:
        data = pd.read_csv(path, sep= '\t', header=None,names=col_names, index_col=None)
    data.dropna(how = 'all', inplace=True, axis=1)
    intensity_adjusted = []
    for index, row in data.iterrows():
        if row['intensity']<0:
            intensity_adjusted.append(0)
        else:
            intensity_adjusted.append(row['intensity'])
    data['intensity_adjusted']=intensity_adjusted
    data = data.query('cs <=9 and cs >= 0.4')
    data.reset_index(inplace=True, drop = True)
    if if_smooth == True:
        data_smoothed = spline_smooth(data['cs'], data['intensity_adjusted'])
        intensity_adjusted = []
        for index, row in data_smoothed.iterrows():
            if row['cs'] >=0.4 and row['cs']<=3.3:
                intensity_adjusted.append(row['intensity_adjusted'])
            elif row['cs']<=4.7 and row['cs']>=3.35:
                intensity_adjusted.append(row['intensity_adjusted'])
            elif row['cs']<=9 and row['cs']>=5:
                intensity_adjusted.append(row['intensity_adjusted'])
            else:
                intensity_adjusted.append(0)
        data_smoothed['intensity_adjusted']=intensity_adjusted
        data_smoothed.sort_values(by = 'cs', ascending=True, inplace=True)
        data_smoothed.reset_index(inplace=True, drop=True)
        return (data_smoothed, data)
    else:

        return(data)
def spline_smooth(x_points, y_points, n = 5):
    tck = interpolate.splrep(x_points, y_points)
    x_fit = np.arange(0.4, 9, 0.000275)
    y_fitted = interpolate.splev(x_fit, tck)
    for i in range(0, len(y_fitted)):
        if y_fitted[i]<0:
            y_fitted[i]=0
    df = pd.DataFrame(zip(x_fit, y_fitted), columns=['cs', 'intensity_adjusted'])
    weights = np.arange(1, n + 1)
    wmas = df['intensity_adjusted'].rolling(n).apply(lambda x: np.dot(x, weights) /
                                                               weights.sum(), raw=True).to_list()
    df['intensity_adjusted'] = wmas
    df.dropna(subset = ['intensity_adjusted'], inplace=True)
    df.reset_index(inplace=True, drop=True)
    return(df)
def nmr_stack(profile1_input, profile2_input, label1 ='fraction',label2='simulated_signal',intensity_col = 'intensity_adjusted',vline_location_upper = None,cs_start =None, cs_end = None, shift2 = None, save_path = None):
    profile1 = profile1_input.copy()
    profile2 = profile2_input.copy()
    fig = plt.figure(figsize = (10, 8))#43
    plt.subplots_adjust()
    ax = fig.add_subplot()
    if shift2 is not None:
        profile2['cs']=profile2['cs']+shift2
    if cs_start is not None and cs_end is not None:
        profile1 = quick_search_values(profile1, 'cs', cs_start, cs_end, ifsorted=False)
        profile2 = quick_search_values(profile2, 'cs', cs_start, cs_end, ifsorted=False)
        # ax.set_xlim(cs_start,cs_end)
    # intensity_normalized1 = profile1[intensity_col] / profile1[intensity_col].max()  * 100.0
    # intensity_normalized2 = -(profile2[intensity_col] / profile2[intensity_col].max()  * 100.0)
    sns.lineplot(x= profile1['cs'], y = profile1[intensity_col], color = 'blue', label = label1)
    sns.lineplot(x= profile2['cs'], y = profile2[intensity_col], color = 'red', label = label2)
    ax.legend()
    # intensity_reference =
    if profile2[intensity_col].max()>profile1[intensity_col].max():
        ylim_max = profile2[intensity_col].max()*1.1

    else:
        ylim_max = profile1[intensity_col].max()*1.1
    # ax.set_ylim(0, +ylim_max)
    # p = peaks[350]
    ax.grid(None)
    ax.set_facecolor('none')
    plt.tight_layout()
    if vline_location_upper is not None:
        for loc in vline_location_upper:
            plt.vlines(x = loc  , ymin = 0, ymax =ylim_max, colors = 'orange')
    if save_path is not None:
        plt.savefig(save_path)
def nmr_head_to_tail(profile1_input, profile2_input, intensity_col = 'intensity_adjusted',vline_location_upper = None, vline_location_lower = None,cs_start =None, cs_end = None, shift2 = None):
    profile1 = profile1_input.copy()
    profile2 = profile2_input.copy()
    fig = plt.figure(figsize = (10, 8))#43
    plt.subplots_adjust()
    ax = fig.add_subplot()
    if shift2 is not None:
        profile2['cs']=profile2['cs']+shift2
    if cs_start is not None and cs_end is not None:
        profile1 = quick_search_values(profile1, 'cs', cs_start, cs_end, ifsorted=False)
        profile2 = quick_search_values(profile2, 'cs', cs_start, cs_end, ifsorted=False)
        # ax.set_xlim(cs_start,cs_end)
    intensity_normalized1 = profile1[intensity_col] / profile1[intensity_col].max()  * 100.0
    intensity_normalized2 = -(profile2[intensity_col] / profile2[intensity_col].max()  * 100.0)
    sns.lineplot(x= profile1['cs'], y = intensity_normalized1, color = 'blue')
    sns.lineplot(x= profile2['cs'], y = intensity_normalized2, color = 'red')
    # intensity_reference =
    ax.set_ylim(-105, +105)
    # p = peaks[350]
    ax.grid(None)
    if vline_location_upper is not None:
        for loc in vline_location_upper:
            plt.vlines(x = loc  , ymin = 0, ymax =100, colors = 'orange')
    if vline_location_lower is not None:
        for loc in vline_location_lower:
            plt.vlines(x = loc  , ymin = -100, ymax =0, colors = 'green')
def nmr_plot(profile, intensity_col = 'intensity_adjusted',vline_location_upper = None, cs_start =None, cs_end = None, peak_stat = None):
    fig = plt.figure(figsize = (10, 8))#43
    plt.subplots_adjust()
    ax = fig.add_subplot()
    if cs_start is not None and cs_end is not None:
        profile = quick_search_values(profile, 'cs', cs_start, cs_end, ifsorted=False)
        # ax.set_xlim(cs_start,cs_end)
    # intensity_mixture_normalized = profile[intensity_col] / profile[intensity_col].max()  * 100.0
    sns.lineplot(x= profile['cs'], y = profile[intensity_col], color = 'blue')
    # intensity_reference =
    # ax.set_ylim(0, profile[intensity_col].max()*1.02)
    # p = peaks[350]
    if vline_location_upper is not None:
        for loc in vline_location_upper:
            plt.vlines(x = loc  , ymin  = profile[intensity_col].min(), ymax =profile[intensity_col].max(), colors = 'orange')
    if peak_stat is not None:
        for index, row in peak_stat.iterrows():
            plt.vlines(x = row['peak_apex']  , ymin = 0, ymax =row['intensity'], colors = 'orange')

def get_peak_statistics(data_smoothed, peaks):
    _smoothed = data_smoothed.copy()
    _smoothed.sort_values(by = 'cs', ascending=True, inplace = True)
    cf_start = []
    cf_end = []
    cf_apex = []
    cf_half_start = []
    cf_half_end = []
    peak_intensity = []
    peak_length_full = []
    peak_length_half = []
    peak_area = []
    for i in range(0, len(peaks)):
        full_length, right_edge, left_edge = get_full_width(x = _smoothed['cs'][peaks[i][0]:peaks[i][2]+1], y= _smoothed['intensity_adjusted'][peaks[i][0]:peaks[i][2]+1], height=0.05)
        fwhf, fwhf_ledge, fwhf_redge = get_full_width(x = _smoothed['cs'][peaks[i][0]:peaks[i][2]+1], y= _smoothed['intensity_adjusted'][peaks[i][0]:peaks[i][2]+1], height=0.5)
        cf_start.append(left_edge)
        cf_end.append(right_edge)
        auc_region = quick_search_values(_smoothed, 'cs',left_edge,right_edge, ifsorted=True)
        peak_area.append(np.trapz(auc_region['intensity_adjusted'], dx = step))

        # center =gaussian_estimator(peaks[i], _smoothed['cs'], _smoothed['intensity_adjusted'])
        # cf_apex.append(center)
        cf_apex.append(_smoothed['cs'][peaks[i][1]])
        peak_length_full.append(full_length)
        peak_length_half.append(fwhf)
        peak_intensity.append(_smoothed['intensity_adjusted'][peaks[i][1]])
    peak_statistics = pd.DataFrame(zip(np.arange(len(peaks)), cf_apex,peak_intensity,peak_area,cf_start, cf_end, peak_length_half, peak_length_full), columns=['peak_id', 'peak_apex', 'intensity', 'area','peak_start','peak_end','peak_length_half', 'peak_length_full'])

    return(peak_statistics)
def filter_peaks(peak_statistics):
    _lower =peak_statistics.query('peak_apex <=3.30 & peak_apex >= 0.4')
    _middle = peak_statistics.query('peak_apex<=4.7 & peak_apex >= 3.35')
    _upper = peak_statistics.query('peak_apex <= 9 & peak_apex >=5')
    peak_statistics = pd.concat([_lower, _middle, _upper], axis=0)
    # return(peak_statistics)
    max_intensity = peak_statistics['intensity'].max()
    peak_statistics_major = peak_statistics[peak_statistics['intensity']>max_intensity * 0.01]
    peak_statistics_major.reset_index(inplace = True, drop = True)
    return (peak_statistics_major)
def get_peaks_nmr(intensity_list: np.ndarray) -> list:
    from scipy.signal import find_peaks
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
        peak_list.append(get_edges_nmr(intensity_list, cur_apex_idx))
    return(peak_list)
def get_full_width(x: np.ndarray, y: np.ndarray, height: float = 0.5) -> float:
    height_half_max = np.max(y) * height
    index_max = np.argmax(y)
    x_low = np.interp(height_half_max, y[:index_max+1], x[:index_max+1])
    x_high = np.interp(height_half_max, np.flip(y[index_max:]), np.flip(x[index_max:]))

    return x_high - x_low, x_high, x_low

def get_edges_nmr(intensity_list, cur_apex_idx, edge_cof = 0.01):
    current_apex_intensity = intensity_list[cur_apex_idx]
    left_edge_int = intensity_list[cur_apex_idx]
    # check if left of the apex is less than current apex
    for i in range(1, cur_apex_idx+1):

        if intensity_list[cur_apex_idx-i]<=left_edge_int and left_edge_int>=edge_cof*current_apex_intensity:
            left_edge_int = intensity_list[cur_apex_idx-i]
            left_edge_idx = cur_apex_idx-i
        else:
            left_edge_idx = cur_apex_idx-i+1
            break
    # going right
    right_edge_int = intensity_list[cur_apex_idx]
    for i in range(1, len(intensity_list)-cur_apex_idx):
        if intensity_list[cur_apex_idx+i]<=right_edge_int and right_edge_int >=edge_cof*current_apex_intensity:
            right_edge_int = intensity_list[cur_apex_idx+i]
            right_edge_idx = cur_apex_idx+i
        else:
            right_edge_idx = cur_apex_idx+i-1
            break
    return((left_edge_idx, cur_apex_idx, right_edge_idx))
def cluster_peaks(peak_major, step = 0.02):
    peak_major_sorted = peak_major.copy()
    peak_major_sorted.sort_values(by = 'intensity', inplace=True, ascending=False)
    cluster_id = 0
    peak_clustered = pd.DataFrame()
    while len(peak_major_sorted)>0:
        current_cluster = quick_search_values(peak_major_sorted,'peak_apex', peak_major_sorted.iloc[0]['peak_apex']-step, peak_major_sorted.iloc[0]['peak_apex']+step, ifsorted=False)
        current_cluster['cluster_id']=cluster_id
        cluster_id+=1
        peak_clustered = pd.concat([peak_clustered, current_cluster])
        peak_major_sorted.drop(list(current_cluster.index.values), inplace=True)
    return peak_clustered
import math
def get_shift_statistics(cluster_statistics, std_profile, fraction_profile):
    step = 0.000275
    peak_cluster_statistics = pd.DataFrame()
    for cluster in cluster_statistics['cluster_id'].unique():
        current_cluster = string_search(cluster_statistics, 'cluster_id', cluster)
        left_edge = current_cluster['peak_start'].min()-0.01
        right_edge = current_cluster['peak_end'].max()+0.01
        std_region = quick_search_values(std_profile, 'cs', left_edge, right_edge, ifsorted=False)
        fraction_region = quick_search_values(fraction_profile, 'cs', left_edge, right_edge, ifsorted=False)
        corr_forward = sm.tsa.stattools.ccf(fraction_region['intensity_adjusted'], std_region['intensity_adjusted'], adjusted=False)

        corr_backward = sm.tsa.stattools.ccf(std_region['intensity_adjusted'],fraction_region['intensity_adjusted'] , adjusted=False)
        if np.max(corr_forward)>np.max(corr_backward):
            lag = np.argmax(corr_forward)
            shift = step*lag
            if math.isinf(np.max(corr_forward))==False:
                max_cor = corr_forward[lag]
            else:
                max_cor = 0
        else:
            lag = np.argmax(corr_backward)
            shift = step*lag*(-1)
            if math.isinf(np.max(corr_backward))==False:
                max_cor = corr_backward[lag]
            else:
                max_cor = 0
        current_cluster['shift']=shift
        current_cluster['corr']=max_cor
        split = len(current_cluster[current_cluster['area']>0.4*current_cluster['area'].max()])
        current_cluster['split']=split

        # if len(current_cluster[current_cluster['area']>0.8*current_cluster['area'].max()])==1:
        #     coe = np.log10(1.1)
        # else:
        #     coe  = np.log10(len(current_cluster[current_cluster['area']>0.8*current_cluster['area'].max()]))
        if split <2:
            coe = 0
        else:
            coe = 1
        # coe = np.log10(split)
        # current_cluster['score']=max_cor*np.log10(coe)*np.log10(current_cluster['area'].median())
        current_cluster['score']=coe*(max_cor*np.log10(current_cluster['area'].max()))
        peak_cluster_statistics = pd.concat([peak_cluster_statistics, current_cluster], ignore_index=True)
    return peak_cluster_statistics
def get_cluster_statistics(profile_smoothed):
    peaks = get_peaks_nmr(profile_smoothed['intensity_adjusted'])
    peak_stats = get_peak_statistics(profile_smoothed, peaks)
    peak_major = filter_peaks(peak_stats)
    peak_clustered = cluster_peaks(peak_major)
    return(peak_clustered)
def get_all_peaks(profile_smoothed):
    peaks = get_peaks_nmr(profile_smoothed['intensity_adjusted'])
    peak_stats = get_peak_statistics(profile_smoothed, peaks)
    peak_stats=peak_stats[peak_stats['area']>0]
    peak_stats.sort_values(by = 'peak_apex', inplace = True, ascending=True)
    return peak_stats
    # current_sum = 0
    # peak_statistics_major = pd.DataFrame()
    # for index, row in peak_statistics_combined.iterrows():
    #     if current_sum <= peak_statistics_combined['intensity'].sum()*coef:
    #         peak_statistics_major = pd.concat([peak_statistics_major, pd.DataFrame([row])], axis = 0)
    #         current_sum = current_sum + row['intensity']
    #     else:
    #         break
    # return(peak_statistics_major)
def get_model_cluster(std_smoothed, fraction_smoothed, if_all = True):
    cluster_statistics = get_cluster_statistics(std_smoothed)
    cluster_shift_statistics = get_shift_statistics(cluster_statistics, std_profile=std_smoothed, fraction_profile=fraction_smoothed)
    cluster_shift_statistics.sort_values(by = 'score', ascending=False, inplace=True)
    model_cluster = cluster_shift_statistics[cluster_shift_statistics['score']==cluster_shift_statistics['score'].max()]
    if if_all ==True:
        return(model_cluster, cluster_shift_statistics)
    else:
        return (model_cluster)
def get_matched_statistics(model_cluster, fraction_smoothed):
    spike_5_peaks_all = get_all_peaks(fraction_smoothed)
    ratio = []
    matched_peaks = pd.DataFrame()
    for index, row in model_cluster.iterrows():
        diff = abs(spike_5_peaks_all['peak_apex']-(row['peak_apex']+row['shift']))

        matched_peak = spike_5_peaks_all.iloc[np.argmin(diff)]
        matched_peaks = pd.concat([matched_peaks, pd.DataFrame([matched_peak])])
        ratio.append(matched_peak['intensity']/row['intensity'])
    matched_peaks.columns = 'fraction_'+matched_peaks.columns
    matched_peaks.reset_index(inplace=True, drop=True)
    model_cluster.reset_index(inplace=True, drop=True)
    matched_statistics = pd.concat([model_cluster, matched_peaks], axis=1)
    matched_statistics.insert(0, 'ratio', ratio)
    return(matched_statistics)

# below are deprecated
# def determine_noise_level(inputt, quantile = 0.1):
#     nmr_data = inputt.copy()
#     x = np.array(inputt['cs'])
#     y = np.array(inputt['intensity_smoothed'])
#     dri = []
#     for i in range(0, len(x)):
#         if i -2<0 or i+2>len(x)-1:
#             dri.append(np.NAN)
#         else:
#             dri_temp = -2*y[i-2]-y[i-1]+y[i+1]+2*y[i+2]
#             dri.append(dri_temp)
#     dri_abs = np.abs(dri)
#     nmr_data['dri']=dri
#     nmr_data['dri_abs']=dri_abs
#     # return(nmr_data)
#     # dri = np.diff(inputt['intensity_smoothed'])/np.diff(inputt['cs'])
#     # dri = np.gradient(y,x)
#     # dri_abs = np.abs(dri)
#     # dri_abs = np.append(np.NAN, dri_abs)
#     nmr_data_good_region = nmr_data[nmr_data['cs_label']!='Removing']
#     nmr_data_good_region.reset_index(inplace = True, drop = True)
#     noise_level = nmr_data_good_region[nmr_data_good_region['dri_abs']>nmr_data_good_region['dri_abs'].quantile(0.99)]['dri_abs'].median()
#     # return(nmr_data_good_region)
#     # nmr_data_noise = nmr_data_good_region[nmr_data_good_region['dri_abs']<nmr_data_good_region['dri_abs'].quantile(quantile)]
#     return(nmr_data_good_region, noise_level)
# def wma(nmr_data, column='intensity_adjusted', n=2):
#
#     intensity_adjusted = []
#     for index, row in nmr_data.iterrows():
#         if row['intensity']<0:
#             intensity_adjusted.append(0)
#         else:
#             intensity_adjusted.append(row['intensity'])
#     nmr_data['intensity_adjusted']=intensity_adjusted
#     weights = np.arange(1, n + 1)
#     wmas = nmr_data[column].rolling(n).apply(lambda x: np.dot(x, weights) /
#                                                  weights.sum(), raw=True).to_list()
#
#     nmr_data['intensity_adjusted'] = wmas
#
#     # nmr_data['intensity_adjusted']=intensity_adjusted
#     nmr_data.dropna(subset = ['intensity_adjusted'], inplace=True)
#     nmr_data.reset_index(inplace=True, drop=True)
#         # df.dropna(subset = ['intensity_smoothed'], inplace=True)
#         # df.reset_index(inplace=True, drop=True)
#     return nmr_data

# def gaussian_estimator(
#         peak: tuple,
#         mz_array: np.ndarray,
#         int_array: np.ndarray
# ) -> float:
#     start, center, end = peak
#
#     m1, m2, m3 = mz_array[center - 1], mz_array[center], mz_array[center + 1]
#     i1, i2, i3 = int_array[center - 1], int_array[center], int_array[center + 1]
#     # print(m1,m2,m3)
#     if i1 == 0:  # Case of sharp flanks
#         m = (m2 * i2 + m3 * i3) / (i2 + i3)
#     elif i3 == 0:
#         m = (m1 * i1 + m2 * i2) / (i1 + i2)
#     else:
#         l1, l2, l3 = np.log(i1), np.log(i2), np.log(i3)
#         m = (
#                 ((l2 - l3) * (m1 ** 2) + (l3 - l1) * (m2 ** 2) + (l1 - l2) * (m3 ** 2))
#                 / ((l2 - l3) * (m1) + (l3 - l1) * (m2) + (l1 - l2) * (m3))
#                 * 1
#                 / 2
#         )
#
#     return m
# def label_region(inputt, cs_col = 'cs'):
#     nmr_data = inputt.copy()
#     labels = []
#     for index, row in nmr_data.iterrows():
#         if row[cs_col]>5 and row[cs_col]<9:
#             labels.append('Upper')
#         elif row[cs_col]>0.4 and row[cs_col]<4.7:
#             labels.append('Low')
#         else:
#             labels.append('Removing')
#     nmr_data['cs_label']=labels
#     return(nmr_data)

#