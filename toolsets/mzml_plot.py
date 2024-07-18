import os
import sys
import numpy as np
import pymzml
from pymzml.plot import Factory
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
def TIC(mzml_file):
    """
    Plots a chromatogram for the given mzML file. File is saved as
    'chromatogram_<mzml_file>.html'.

    usage:

        ./plot_chromatogram.py <path_to_mzml_file>

    """
    run = pymzml.run.Reader(mzml_file)
    mzml_basename = os.path.basename(mzml_file)
    pf = Factory()
    pf.new_plot()
    pf.add(run["TIC"].peaks(), color=(0, 0, 0), style="lines", title=mzml_basename)
    pf.save(
        "chromatogram_{0}.html".format(mzml_basename),
        layout={"xaxis": {"title": "Retention time"}, "yaxis": {"title": "TIC"}},
    )
    return
def EIC(mzml_file, ms2_value, rt_start = -1, rt_end = -1, vlines=False, vlines_location = []):

    run = pymzml.run.Reader(mzml_file,

        MS_precisions ={
            1 : 10e-6,
            2 : 20e-6
        }
        )
    time_dependent_intensities = []

    # ms2_value = 271.059450
    for spectrum in tqdm(run, total=run.get_spectrum_count()):
        if spectrum.ms_level == 1:
            has_peak_matches = spectrum.has_peak(ms2_value)
            if has_peak_matches != []:
                for mz, I in has_peak_matches:
                    time_dependent_intensities.append(
                        [spectrum.scan_time_in_minutes(), I, mz]
                    )
    rt = []
    intensity = []
    for scan in time_dependent_intensities:
        # break
        rt.append(scan[0])
        intensity.append(scan[1])

    fig, ax = plt.subplots(
    figsize = (9, 5)
                      )
    try:
        ax= sns.lineplot(x=rt, y = intensity)
    except:
        return(rt, intensity)


    if rt_start != -1 and rt_end != -1:
        ax.set_xlim(rt_start, rt_end)
        ax.set_ylim(0, np.max(intensity)+100)
    if vlines == True and vlines_location != []:
        for position in vlines_location:
            plt.axvline(x = position, color = 'red', label = row['ms1_precursor_intensity'])
        

    ax.grid(False)
    # return()
    

