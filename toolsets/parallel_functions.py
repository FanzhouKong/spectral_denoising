import pandas as pd
import numpy as np
import toolsets.spectra_operations as so
from matplotlib import rcParams
# from mimas.external.features_by_alphapept.load_mzml_data import load_mzml_data
# import toolsets.spectra_operations as so
from toolsets.spectra_operations import break_spectra, pack_spectra
import matplotlib.pyplot as plt

import seaborn as sns
from toolsets.spectra_plotter import chop_msms, _extract_MS1
import plotly.express as px
import pymzml
def plus(a,b):
    return(a+b)
def _extract_ms1_intensity_p(ms1, idx, pmz):
   mass_temp, intensity_temp = so.break_spectra(ms1.iloc[idx]['peaks'])

   index_start = np.searchsorted(mass_temp, pmz-0.005,side = 'left')
   index_end = np.searchsorted(mass_temp, pmz+0.005,side = 'right')
   return(ms1.iloc[idx]['rt'],np.sum(intensity_temp[index_start: index_end]))
def _extract_mzml_info(index, spec_list) -> tuple:
    """Extract basic MS coordinate arrays from a dictionary.

    Args:
        input_dict (dict): A dictionary obtained by iterating over a Pyteomics mzml.read function.

    Returns:
        tuple: The rt, masses, intensities, ms_order, prec_mass, mono_mz, charge arrays retrieved from the input_dict.
            If the `ms level` in the input dict does not equal 2, the charge, mono_mz and prec_mass will be equal to 0.

    """
    input_dict = spec_list[index]
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
    prec_mass = mono_mz * abs(charge) - charge * 1.00727646687
    return prec_mass
