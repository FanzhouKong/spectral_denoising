import numpy as np
import pymzml
import pandas as pd
import toolsets.spectra_operations as so
from toolsets.search import string_search
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

def mzml_preprocess(mzml, if_mix = False):
    if if_mix == True:
        mzml_data = load_mzml_data(mzml+'.mzml')
    else:
        mzml_data = load_mzml_data(mzml)
    ms2 = string_search(mzml_data, 'ms_level', 2)
    ms2.reset_index(inplace=True, drop=True)
    ms2.sort_values(by='precursor_mz', inplace=True)
    ms1 = string_search(mzml_data, 'ms_level', 1)
    ms1.reset_index(inplace=True, drop=True)
    ms1_intensity = []
    peak_purity = []
    for index, row in ms2.iterrows():
        ms1_scan = string_search(ms1, 'cycle', row['cycle']).iloc[0]
        precursor_intensity = _extract_ms1_intensity(ms1_scan['peaks'], row['precursor_mz']-0.005, row['precursor_mz']+0.005)
        isolation_window_intensity = _extract_ms1_intensity(ms1_scan['peaks'], row['isolation_window'][0], row['isolation_window'][1])
        ms1_intensity.append(precursor_intensity)
        peak_purity.append(precursor_intensity/isolation_window_intensity)
    ms2['ms1_intensity'] =ms1_intensity
    ms2['peak_purity']=peak_purity
    return(ms2)
def _extract_ms1_intensity(peaks, mz_lower, mz_upper):
    mass_temp, intensity_temp = so.break_spectra(peaks)
    index_start = np.searchsorted(mass_temp, mz_lower,side = 'left')
    index_end = np.searchsorted(mass_temp, mz_upper,side = 'right')
    return(np.sum(intensity_temp[index_start: index_end]))





def load_mzml_data(file: str, n_most_abundant=500, nested_array = False) -> tuple:
    """Load data from an mzml file as a dictionary.

    Args:
        filename (str): The name of a .mzml file.
        n_most_abundant (int): The maximum number of peaks to retain per MS2 spectrum.
        callback (callable): A function that accepts a float between 0 and 1 as progress. Defaults to None.

    Returns:
        tuple: A dictionary with all the raw data, a string with the acquisition_date_time and a string with the vendor.

    """
    import pymzml
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

    scan_list = []
    rt_list = []
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
    for spec in pymzml.run.Reader(file, obo_version="4.1.33"):
        rt, masses, intensities, ms_order, prec_mass, mono_mz, charge, (prec_windows_lower, prec_windows_upper) = _extract_mzml_info(spec)

        # Remove zero intensities
        to_keep = intensities > 0
        masses = masses[to_keep]
        intensities = intensities[to_keep]
        if ms_order == 1:
            cycle = cycle+1
        cycles.append(cycle)
        # Sort by m/z
        if not np.all(masses[:-1] <= masses[1:]):
            order = np.argsort(masses)
            masses = masses[order]
            intensities = intensities[order]

        # Only keep top n_most_abundant peaks
        if ms_order == 2 and len(masses) > n_most_abundant:
            sortindex = np.argsort(intensities)[::-1][:n_most_abundant]
            sortindex.sort()
            masses, intensities = masses[sortindex], intensities[sortindex]
        id = id+1

        # mass_list.append(masses)
        # int_list.append(intensities)
        peak_list.append(so.pack_spectra(masses, intensities)) 
        ms_list.append(ms_order)
        prec_mzs_list.append(prec_mass)
        mono_mzs_list.append(mono_mz)
        charge_list.append(charge)
        select_windows_list.append((prec_windows_lower, prec_windows_upper))
        scan_idx.append(cur_idx)
        rt_list.append(rt)
        cur_idx = cur_idx+1


    return pd.DataFrame(
        list(zip(scan_idx, cycles, ms_list,mono_mzs_list, rt_list,charge_list,peak_list,select_windows_list)), 
        columns=['scan_id', 'cycle','ms_level', 'precursor_mz', 'rt','charge', 'peaks', 'isolation_window']
        )

    # scan_list_ms1 = [scan_list[i] for i, _ in enumerate(ms_list) if _ == 1]
    # rt_list_ms1 = [rt_list[i] for i, _ in enumerate(ms_list) if _ == 1]
    # mass_list_ms1 = [mass_list[i] for i, _ in enumerate(ms_list) if _ == 1]
    # int_list_ms1 = [int_list[i] for i, _ in enumerate(ms_list) if _ == 1]
    # ms_list_ms1 = [ms_list[i] for i, _ in enumerate(ms_list) if _ == 1]

    # scan_list_ms2 = [scan_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # rt_list_ms2 = [rt_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # mass_list_ms2 = [mass_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # int_list_ms2 = [int_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # ms_list_ms2 = [ms_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # prec_mass_list2 = [prec_mzs_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # mono_mzs2 = [mono_mzs_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # charge_ms2 = [charge_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # select_windows_ms2 = [select_windows_list[i] for i, _ in enumerate(ms_list) if _ == 2]
    # if nested_array== True:
    #     return(mass_list_ms1, int_list_ms1)
    # prec_mass_list2 = [_calculate_mass(mono_mzs_list[i], charge_list[i])
    #                    for i, _ in enumerate(ms_list) if _ == 2]

    # ms_run = {}
    # # self.clear()
    # ms_run["scan_list_ms1"] = np.array(scan_list_ms1)
    # ms_run["rt_list_ms1"] = np.array(rt_list_ms1)
    # ms_run["mass_list_ms1"] = np.array(mass_list_ms1)
    # ms_run["int_list_ms1"] = np.array(int_list_ms1)
    # ms_run["ms_list_ms1"] = np.array(ms_list_ms1)

    # ms_run["scan_list_ms2"] = np.array(scan_list_ms2)
    # ms_run["rt_list_ms2"] = np.array(rt_list_ms2)
    # ms_run["mass_list_ms2"] = mass_list_ms2
    # ms_run["int_list_ms2"] = int_list_ms2
    # ms_run["ms_list_ms2"] = np.array(ms_list_ms2)
    # ms_run["prec_mass_list2"] = np.array(prec_mass_list2)
    # ms_run["mono_mzs2"] = np.array(mono_mzs2)
    # ms_run["charge2"] = np.array(charge_ms2)
    # ms_run["select_windows_ms2"] = np.array(select_windows_ms2, dtype=np.float32)

    # # Reorganize data: flatten lists
    # ms_run["indices_ms1"] = _index_ragged_list(ms_run["mass_list_ms1"])
    # ms_run["mass_list_ms1"] = np.concatenate(ms_run["mass_list_ms1"])
    # ms_run["int_list_ms1"] = np.concatenate(ms_run["int_list_ms1"])

    # return ms_run
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
def _calculate_mass(mono_mz: float, charge: int) -> float:
    """Calculate the precursor mass from mono mz and charge.

    Args:
        mono_mz (float): mono m/z.
        charge (int): charge.

    Returns:
        float: precursor mass.
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
