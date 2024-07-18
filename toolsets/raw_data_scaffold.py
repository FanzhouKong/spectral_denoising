import pymzml
import numpy as np
import pandas as pd
from toolsets.search import string_search, quick_search_sorted, quick_search_values
import toolsets.helpers as helpers
from toolsets.spectra_operations import break_spectra, pack_spectra
from matplotlib import rcParams
import seaborn as sns
import os
import toolsets.T_rex as trx
import matplotlib.pyplot as plt

def auto_EIC(pmz, ms1, mass_tolerance = 0.005, show_eic = False):
    mass_sorted, intensity_sorted, index_sorted, rt_list = trx.build_index(ms1)
    intensity_list = trx.flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted, mass_error=mass_tolerance)
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
    ax= sns.lineplot(x = rt_list, y = intensity_list, label = 'EIC', color = 'blue')
    max = np.max(intensity_list)
    if base_line_series is not None:
        ax = sns.lineplot(x = rt_list, y = base_line_series, color = 'orange', label = 'baseline', linestyle='--')
        if np.max(base_line_series)>max:
            max = np.max(base_line_series)

    ax.set_ylim(0, max*1.1)
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
            plt.axvline(x = position, color = 'red',linestyle='--')
    if len(vlines_location_2)>0:
        for position in vlines_location_2:
            plt.axvline(x = position, color = 'green',linestyle='--')
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
def read_mzml(mzml_path, parent_dir =  None, rt_max = None,if_mix = False, sp_removal = False):

    if if_mix == True and parent_dir != None:
        mix = mzml_path
        mzml_base_name = helpers.find_files(parent_dir, mix)
        mzml_path = os.path.join(parent_dir, mzml_base_name)
    if if_mix == False and parent_dir is not None:
        mzml_path = os.path.join(parent_dir,mzml_path)
    if mzml_path[-5:]!='.mzML' and mzml_path[-5:]!='.mzml':
        mzml_path = mzml_path+'.mzML'
    ms1_2 = _load_mzml_data(mzml_path, rt_max = rt_max)
    ms2 = string_search(ms1_2, 'ms_level', 2)
    ms1 = string_search(ms1_2, 'ms_level', 1)

    ms1_intensity = np.zeros(len(ms2))
    ms1_rt = np.zeros(len(ms2))
    peak_purity = np.zeros(len(ms2))
    ms1_pmz=np.zeros(len(ms2))
    mz_offsets = np.zeros(len(ms2))
    for index, row in ms2.iterrows():
        ms1_scan = string_search(ms1, 'cycle', row['cycle']).iloc[0]
        pmz_ms1, precursor_intensity, mz_offset= _extract_precursor_intensity(ms1_scan['peaks'], row['precursor_mz'])#search for mass that most close to precursor_mz
        ms1_rt[index]=ms1_scan['rt']
        isolation_window_intensity = _extract_ms1_intensity(ms1_scan['peaks'], row['isolation_window'][0], row['isolation_window'][1] )
        if isolation_window_intensity!= 0:
            peak_purity[index]=(precursor_intensity/isolation_window_intensity)
        else:
            peak_purity[index]=0
        ms1_intensity[index]=precursor_intensity

        ms1_pmz[index]=pmz_ms1
        mz_offsets[index]=mz_offset
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
        ms1['base_name']=mzml_base_name
        ms1['mix']=mix
    else:
        ms2['mix']= os.path.basename(mzml_path[:-5])
        ms1['mix']=os.path.basename(mzml_path[:-5])
        ms2['base_name']=os.path.basename(mzml_path)
        ms1['base_name']=os.path.basename(mzml_path)
    # ms2.sort_values(by = ['ms1_pmz', 'ms1_precursor_intensity'], inplace = True)
    # ms2 =ms2[ms2['peak_purity'] !=0]
    # ms2 = ms2[ms2['mz_offset']<0.5]
    ms2.reset_index(inplace=True, drop=True)

    if sp_removal == True:
        peaks_binned =[]
        for index, row in (ms1.iterrows()):
            peaks_binned.append(np.array(_bin_ms1(row['peaks'],resolution = 60000)).T)
        ms1['peaks_binned']=peaks_binned
    return(ms1, ms2)
def _load_mzml_data(file: str, n_most_abundant=400, rt_max = None) -> tuple:
    """Load data from an mzml file as a dictionary.

    Args:
        filename (str): The name of a .mzml file.
        n_most_abundant (int): The maximum number of peaks to retain per MS2 spectrum.
        callback (callable): A function that accepts a float between 0 and 1 as progress. Defaults to None.

    Returns:
        a pandas dataframe with all the raw data

    """

    id = 1
    all_specs = pymzml.run.Reader(file, obo_version="4.1.33")
    total_spec_count = all_specs.info['spectrum_count']
    # ms_list = np.zeros(total_spec_count,dtype=int)
    ms_list = [None]*total_spec_count
    # mono_mzs_list = np.zeros(total_spec_count)
    mono_mzs_list = [None]*total_spec_count
    mono_mzs_intensity_list = np.zeros((total_spec_count))
    # pmz_intensity_list = np.zeros(total_spec_count)
    pmz_intensity_list = [None]*total_spec_count
    polarity_list = [None]*total_spec_count
    select_windows_list = [None]*total_spec_count
    # cycles = np.zeros(total_spec_count,dtype=int)
    cycles = [None]*total_spec_count
    cycle = -1
    # scan_idx = np.zeros(total_spec_count,dtype=int)
    scan_idx = [None]*total_spec_count
    cur_idx = 0
    peak_list = [None]*total_spec_count
    # rt_list = np.zeros(total_spec_count)
    rt_list = [None]*total_spec_count
    for spec in pymzml.run.Reader(file, obo_version="4.1.33", build_index_from_scratch=True):
        try:
            rt, masses, intensities, ms_order,mono_mz,mono_mz_intensity,polarity, (prec_windows_lower, prec_windows_upper)=_extract_mzml_info(spec)
            if rt_max!= None and rt>rt_max:
                break

            to_keep = intensities > 0
            # return_mass = masses
            try:
                masses = masses[to_keep]
                intensities = intensities[to_keep]
            except:
                continue

            if ms_order == 1:
                cycle = cycle+1
            cycles[cur_idx]=cycle
            if not np.all(masses[:-1] <= masses[1:]):
                order = np.argsort(masses)
                masses = masses[order]
                intensities = intensities[order]

            # Only keep top n_most_abundant peaks
            # if ms_order == 2 and len(masses) > n_most_abundant:
            #     sortindex = np.argsort(intensities)[::-1][:n_most_abundant]
            #     sortindex.sort()
            #     masses, intensities = masses[sortindex], intensities[sortindex]
            peak = np.array([masses, intensities], dtype=np.float64).T
            # peak = pack_spectra(masses, intensities)
            peak_list[cur_idx]=peak
            ms_list[cur_idx]=ms_order
            mono_mzs_list[cur_idx]=mono_mz
            mono_mzs_intensity_list[cur_idx]=mono_mz_intensity
            # pmz_intensity_list[cur_idx]=precursor_intensity
            polarity_list[cur_idx]=polarity
            select_windows_list[cur_idx]=(prec_windows_lower, prec_windows_upper)
            scan_idx[cur_idx]=cur_idx
            rt_list[cur_idx]=rt
            cur_idx = cur_idx+1
        except:
            pass
    return_df = pd.DataFrame({
            'scan_idx':scan_idx,
            'cycle':cycles,
            'ms_level':ms_list,
            'precursor_mz':mono_mzs_list,
            'precursor_intensity':mono_mzs_intensity_list,
            # 'precursor_intensity':pmz_intensity_list,
            'polarity':polarity_list,
            'rt':rt_list,
            'peaks':peak_list,
            'isolation_window':select_windows_list
        })
    return(return_df)






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
    mono_mz = 0
    mono_mz_intensity = 0
    polarity = np.NAN
    # prec_mass = mono_mz = charge = 0
    # if ms_order == 2 and len(input_dict.selected_precursors) > 0:
    if ms_order == 1 and input_dict['isolation window target m/z'] is None:
        polarity = '+'
        if input_dict['negative scan'] is not None:
            polarity = '-'
        elif input_dict['positive scan'] is None:
            raise Exception("Can't determine polarity")
    if ms_order == 2 and input_dict['isolation window target m/z'] is not None:
        polarity = '+'
        if input_dict['negative scan'] is not None:
            polarity = '-'
        elif input_dict['positive scan'] is None:
            raise Exception("Can't determine polarity")
        mono_mz = input_dict['isolation window target m/z']
        if 'i' in input_dict.selected_precursors[0].keys():
            mono_mz_intensity = input_dict.selected_precursors[0]['i']
        else:
            mono_mz_intensity = np.NAN
        # precursor_intensity =input_dict.selected_precursors[0]["i"]
        # prec_mass = _calculate_mass(mono_mz, charge)

        prec_windows_center = mono_mz
        try:
            prec_windows_lower = prec_windows_center-(input_dict["isolation window lower offset"])
            prec_windows_upper = prec_windows_center+(input_dict["isolation window upper offset"])
        except:
            prec_windows_lower = prec_windows_center-0.5
            prec_windows_upper = prec_windows_center+0.5
    else:
        prec_windows_lower, prec_windows_upper = 0., 0.

    return (rt, masses, intensities, ms_order,
            # prec_mass,
            mono_mz,mono_mz_intensity, polarity, (prec_windows_lower, prec_windows_upper))
def auto_TIC(ms1):
    intensity_list = []
    rt_list = []
    for index, row in ms1.iterrows():
        intensity_list.append(np.sum(row['peaks'].T[1]))
        rt_list.append(row['rt'])
    return(rt_list, intensity_list)
def _extract_precursor_intensity(peaks, pmz):
    mass_temp, intensity_temp =peaks.T
    # mass_temp, intensity_temp = so.break_spectra(peaks)
    search_array=np.array(mass_temp)
    index_start, index_end = search_array.searchsorted([pmz-0.5, pmz+0.5])
    mass_temp = mass_temp[index_start:index_end]
    intensity_temp = intensity_temp[index_start:index_end]
    if len(mass_temp)==0:
        return(pmz, 0, np.NAN)
    # ms1_pmz_idx = np.argmax(intensity_temp)
    offsets = [abs(x-pmz) for x in mass_temp]
    ms1_pmz_idx = np.argmin(offsets)
    ms1_pmz = mass_temp[ms1_pmz_idx]
    ms1_pmz_intensity = intensity_temp[ms1_pmz_idx]
    mz_offset = abs(ms1_pmz-pmz)
    return(ms1_pmz, ms1_pmz_intensity, mz_offset)
def _extract_ms1_intensity(peaks, mz_lower, mz_upper):
    mass_temp, intensity_temp =peaks.T
    search_array=np.array(mass_temp)
    index_start, index_end = search_array.searchsorted([mz_lower, mz_upper])
    return(np.sum(intensity_temp[index_start: index_end]))
def _bin_ms1(peaks, resolution = 60000):
    mass, intensity = peaks.T
    mass_binned = np.zeros(len(mass))
    intensity_binned = np.zeros(len(mass))
    counter = 0
    while len(mass)>0:
        seed_idx = np.argmax(intensity)
        seed_mass = mass[seed_idx]
        res = seed_mass/resolution
        idx_left, idx_right = mass.searchsorted([seed_mass-3*res, seed_mass+3*res])
        mass_binned[counter]=seed_mass
        intensity_binned[counter]=intensity[seed_idx]
        mass = np.concatenate((mass[0:idx_left], mass[idx_right:]))
        intensity = np.concatenate((intensity[0:idx_left], intensity[idx_right:]))
        counter = counter+1#make sure this always the last line!
    mass_binned = mass_binned[0:counter]
    intensity_binned = intensity_binned[0:counter]
    order = np.argsort(mass_binned)
    mass_binned = mass_binned[order]
    intensity_binned = intensity_binned[order]

    return((mass_binned, intensity_binned))