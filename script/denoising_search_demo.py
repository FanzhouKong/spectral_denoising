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
    result = sd.denoising_search(query_peak, query_pmz, reference_library)
    entropy_similarity = result.iloc[0]['entropy_similarity']
    denoised_similarity = result.iloc[0]['denoised_similarity']
    print(f'the entropy similarity is {entropy_similarity} and the denoised similarity is {denoised_similarity} ')
    print('start denoising search in batch mode')
    results = sd.denoising_search_batch(query_spectra['peaks'], query_spectra['precursor_mz'],reference_library)

    result = results[6]
    annotation = result.iloc[0]['name']
    denoised_similarity = result.iloc[0]['denoised_similarity']
    print(f'the denoising search annotation is {annotation}, with the denoised similarity is {denoised_similarity}')
if __name__ == "__main__":
    main()
    print('passed tester!')