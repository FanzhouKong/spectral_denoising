import numpy as np
import sys
import os
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt

import spectral_denoising as sd
from spectral_denoising.spectral_operations import *

sys.path.insert(0, '..')

def main():
    query_data = sd.read_msp('../sample_data/noisy_spectra.msp').iloc[0]
    reference_data= sd.read_msp('../sample_data/clean_spectra.msp').iloc[0]
    query_spectrum, query_smiles, query_adduct, query_pmz = query_data['peaks'], query_data['smiles'], query_data['adduct'], query_data['precursor_mz']
    reference_spectrum = reference_data['peaks']
    denoised_spectrum = sd.spectral_denoising(query_spectrum, query_smiles, query_adduct) 
    print(f'the denoised spectral has similarity of {entropy_similairty(denoised_spectrum, reference_spectrum, pmz = query_pmz):.2f}')
    print('start spectral denoising in batch mode')
    query_data = sd.read_msp('../sample_data/noisy_spectra.msp')
    reference_data= sd.read_msp('../sample_data/clean_spectra.msp').iloc[0]

    pmz, reference_spectra = reference_data['precursor_mz'],reference_data['peaks']
    query_spectra,query_smiles,query_adducts = query_data['peaks'],query_data['smiles'],query_data['adduct']

    denoised_spectra = sd.spectral_denoising_batch(query_spectra,query_smiles,query_adducts)
    print(f'the frist denoised spectra has similarity of {entropy_similairty(denoised_spectra[0], reference_spectrum, pmz = pmz):.2f}')
if __name__ == "__main__":
    main()
    print('passed tester!')


