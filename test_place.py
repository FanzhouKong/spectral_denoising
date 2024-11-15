import numpy as np
import spectral_denoising as sd
from spectral_denoising.spectral_operations import *
from spectral_denoising.chem_utils import *

smiles = 'O=c1nc[nH]c2nc[nH]c12'
formula = 'C5H4N4O'
adduct = '[M+Na]+'
pmz = calculate_precursormz(adduct,smiles)
peak = np.array([
                [48.992496490478516 ,154.0],
                 [63.006099700927734, 265.0],
                 [63.99686813354492, 663.0],
                 [65.9947509765625, 596.0],
                 [79.02062225341797, 521.0],
                 [81.01649475097656, 659.0],
                 ], dtype = np.float32)
print(f'the spectrum entropy of raw spectrum is {spectral_entropy(peak):.2f}, the normalized entropy of raw spectrum is {normalized_entropy(peak):.2f}')
# alternatively, you can store mass and intensity in separate arrays, and use pack_spectrum(mass, intensity) to get the peaks array
# e.g.mass,intensity = [48.992496490478516, 63.006099700927734, 79.02062225341797], [154.0, 265.0, 521.0]
# peak = pack_spectrum(mass, intensity)

# generate some noise ions and add it to the peaks
from spectral_denoising.noise import *
peak_with_noise= sd.read_msp('sample_data/noisy_spectra.msp').iloc[0]['peaks']
peak_denoised = sd.spectral_denoising(peak_with_noise, smiles, adduct)
# peak_denoised = sd.spectral_denoising(peak_with_noise, formula, adduct) # all function would also work if you choose to use formula information instead of smiles information
# use head_to_tail_plot to visualize the spectra, only in jupyter notebook
# sd.head_to_tail_plot(peaks_with_noise,peaks ,pmz)
print(f'the spectrum entropy of contaminated spectrum is {spectral_entropy(peak_with_noise):.2f}, the normalized entropy of contaminated spectrum is {normalized_entropy(peak_with_noise):.2f}')
print(f'the entropy similarity of contaminated spectrum and the raw spectrum is {entropy_similairty(peak_with_noise,peak,  pmz = pmz):.2f}')

# perform spectral denosing and compare against the raw spectrum
peak_denoised = sd.spectral_denoising(peak_with_noise, smiles, adduct)
# peak_denoised = sd.spectral_denoising(peak_with_noise, formula, adduct) # this two would both work, if you wish to use formula information instead of SMILES information
print(f'the entropy similarity of denoised spectrum and the raw spectrum is {entropy_similairty(peak_denoised, peak, pmz = pmz):.2f}')
# use head_to_tail_plot to visualize the spectra, only in jupyter notebook
# sd.head_to_tail_plot(peaks_denoised,peaks ,pmz)