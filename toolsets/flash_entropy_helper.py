from ms_entropy import FlashEntropySearch
import toolsets.spectra_operations as so
import numpy as np
import pprint
spectral_library = [{
    "id": "Demo spectrum 1",
    "precursor_mz": 150.0,
    "peaks": [[100.0, 1.0], [101.0, 1.0], [103.0, 1.0]]
    }]
def flash_entropy_simple(peak1,peak2, pmz):
    if so.num_peaks(peak1)==0 or so.num_peaks(peak2)==0:
        # print('i am in if')
        return(0)
    spectral_library = [{
        "id": "peak2",
        "precursor_mz": pmz,
        "peaks": format_peaks(peak2)
    }]
    noise_level = 0
    query_spectrum = {"precursor_mz": pmz,
                      "peaks": np.array(format_peaks(peak1), dtype=np.float32)}
    entropy_search = FlashEntropySearch(max_ms2_tolerance_in_da=0.02)
    spectral_library = entropy_search.build_index(spectral_library)
    entropy_similarity = entropy_search.search(
        precursor_mz=query_spectrum['precursor_mz'], peaks=query_spectrum['peaks'], noise_threshold=noise_level, min_ms2_difference_in_da = 0.02)
    return(entropy_similarity['identity_search'][0])
def flash_entropy(peak1, pmz1, peak2, pmz2):
    spectral_library = [{
    "id": "peak2",
    "precursor_mz": pmz2,
    "peaks": format_peaks(peak2)
    }]
    query_spectrum = {"precursor_mz": pmz1,
                  "peaks": np.array(format_peaks(peak1), dtype=np.float32)}
    entropy_search = FlashEntropySearch()
    spectral_library = entropy_search.build_index(spectral_library)
    entropy_similarity = entropy_search.search(
    precursor_mz=query_spectrum['precursor_mz'], peaks=query_spectrum['peaks'])
    return(entropy_similarity)
def format_peaks(peaks_string):

    mass, intensity = so.break_spectra(peaks_string)
    lst = []
    for i in range(len(mass)):
        lst.append([mass[i], intensity[i]])
    return(lst)




entropy_search = FlashEntropySearch()
entropy_search.build_index(spectral_library)