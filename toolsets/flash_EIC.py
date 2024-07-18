import numpy as np
import pymzml
import pandas as pd
import itertools
def read_mzml(mzml_path, rt_max = None):

    ms1_2 = _load_mzml_data(mzml_path, rt_max = rt_max)
    ms1 = string_search(ms1_2, 'ms_level',1)
    ms2 = string_search(ms1_2, 'ms_level',2)

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
    ms2['ms1_pmz']=ms1_pmz
    ms2['ms1_rt']=ms1_rt
    ms2['ms1_precursor_intensity'] =ms1_intensity
    ms2['peak_purity']=peak_purity
    ms2['mz_offset']=mz_offsets
    ms2.reset_index(inplace=True, drop=True)
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
            rt, masses, intensities, ms_order,mono_mz,polarity, (prec_windows_lower, prec_windows_upper)=_extract_mzml_info(spec)
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
        # 'precursor_intensity':pmz_intensity_list,
        'polarity':polarity_list,
        'rt':rt_list,
        'peaks':peak_list,
        'isolation_window':select_windows_list
    })
    return(return_df)



def string_search(data, column_name,item, reset_index = True,reverse = False):
    if reverse == False:
        _data= data[data[column_name].to_numpy() == item]
        # return data[data[column_name].to_numpy() == item]
    else:
        _data= data[data[column_name].to_numpy() != item]
    if reset_index == True:
        _data.reset_index(inplace= True, drop = True)
    return(_data)
        # return data[data[column_name].to_numpy() != item]



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
        # precursor_intensity =input_dict.selected_precursors[0]["i"]
        # prec_mass = _calculate_mass(mono_mz, charge)
        try:
            prec_windows_center = (input_dict['isolation window target m/z'])
        except:
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
            mono_mz,polarity, (prec_windows_lower, prec_windows_upper))

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
def flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted, mass_error=0.005):
    index_start, index_end = mass_sorted.searchsorted([pmz-mass_error, pmz+mass_error])
    index_range = index_sorted[index_start:index_end]

    intensity_range = intensity_sorted[index_start:index_end]

    intensity_list = np.zeros(np.max(index_sorted)+1)
    for idx in range(0,len(index_range)):
        intensity_list[index_range[idx]]= intensity_list[index_range[idx]]+intensity_range[idx]
    return(intensity_list)


if __name__ == "__main__":
    mzml_path = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/benchmarking_dataset/mzml/SA1.mzML'
    ms1, ms2 = read_mzml(mzml_path)
    mass_sorted, intensity_sorted, index_sorted, rt_list = build_index(ms1)
    pmz= 371.112
    ion_trace = flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted, mass_error=0.005)
    print(ion_trace[0:10])