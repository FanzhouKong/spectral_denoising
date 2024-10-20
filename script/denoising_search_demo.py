import warnings
warnings.filterwarnings("ignore")

import spectral_denoising as sd
import pandas as pd
import sys
import os
sys.path.insert(0, '..')
def main():
    query_spectra= sd.read_msp('../sample_data/query_spectra.msp')
    reference_library =sd.read_msp('../sample_data/reference_library.msp')
    query_peak, query_pmz = query_spectra.iloc[0]['peaks'], query_spectra.iloc[0]['precursor_mz']
    print(sd.denoising_search(query_peak, query_pmz, reference_library))
    print('start denoising search in batch mode')
    results = sd.denoising_search_batch(query_spectra['peaks'], query_spectra['precursor_mz'],reference_library)
    print(results[6])
if __name__ == "__main__":
    main()