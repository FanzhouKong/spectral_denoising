# Spectral denoising and denoising search

Test the package at: https://pypi.org/project/spectral-denoising/


If you have any questions, feel free to send me E-mails: fzkong@ucdavis.edu. If you find this package useful, please consider citing the following papers:

> . **Denoising Search doubles the number of metabolite and exposome annotations in human plasma using an Orbitrap Astral mass spectrometer**, Res Sq, 2024). [https://doi.org/10.1038/s41592-023-02012-9]

# Note on 2024-10-19

**Please use Python version between 3.8 to 3.12 for this package to work. RDkit currently does not have a distribution compitable to python 3.13!!!!!**

# Project information

Remove noise ions from MS/MS spectra has been tackled for years by mass spectrometrists. Noise ions in MS/MS spectra are largely categorized as 1. electronic noises and 2. chemical noises. 
In this project, we aim to eliminate both chemical noise and electronic noises for improving high-confidence compound identification. 
Integrating such process into spectra matching process, we developed denoising search, which psudo-denoise spectra based on molecular information fetched from reference databases.
This project also provides useful tools to read, write, visualize and compare spectra.

# How to use this package

This repository in Python. A python version >= 3.8 is preferred. 

Detailed documentation can be found at: https://spectral-denoising.readthedocs.io/en/latest/index.html

### Installation

```bash
pip install spectral-denoising
```

### Usage of Classic spectral denoising (electronic denoising and chemical denoising)
The demo data used here can be found under sample_data directory.
#### Simple usage on single spectra
```python
import numpy as np
import spectral_denoising as sd
from spectral_denoising.spectral_operations import *
from spectral_denoising.chem_utils import *

smiles = 'O=c1nc[nH]c2nc[nH]c12'
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
# use head_to_tail_plot to visualize the spectra, only in jupyter notebook
# sd.head_to_tail_plot(peaks_with_noise,peaks ,pmz)
print(f'the spectrum entropy of contaminated spectrum is {spectral_entropy(peak_with_noise):.2f}, the normalized entropy of contaminated spectrum is {normalized_entropy(peak_with_noise):.2f}')
print(f'the entropy similarity of contaminated spectrum and the raw spectrum is {entropy_similairty(peak_with_noise,peak,  pmz = pmz):.2f}')

# perform spectral denosing and compare against the raw spectrum
peak_denoised = sd.spectral_denoising(peak_with_noise, smiles, adduct)
print(f'the entropy similarity of denoised spectrum and the raw spectrum is {entropy_similairty(peak_denoised, peak, pmz = pmz):.2f}')
# use head_to_tail_plot to visualize the spectra, only in jupyter notebook
# sd.head_to_tail_plot(peaks_denoised,peaks ,pmz)
```
#### Spectral denoising on the all spectra from .msp file
```python
import spectral_denoising as sd
query_data = sd.read_msp('sample_data/noisy_spectra.msp')
query_peaks,query_smiles,query_adduct, query_pmz = query_data['peaks'],query_data['smiles'],query_data['adduct'], query_data['precursor_mz'] 
desnoied_peaks = sd.spectra_denoising_batch(query_peaks,query_smiles,query_adduct) # this will return all denoised spectra in a list
```

### Usage of Denoising search
#### Denoising search on a single spectrum against reference library
```python
import spectral_denoising as sd
query_spectra= sd.read_msp('sample_data/query_spectra.msp')
reference_library =sd.read_msp('sample_data/reference_library.msp')
query_spectrum, query_pmz = query_spectra.iloc[0]['peaks'], query_spectra.iloc[0]['precursor_mz'] # just the first spectrum
result = sd.denoising_search(query_spectrum, query_pmz, reference_library)
# result will return all precursor candidates of the query spectrum, each with entropy similarities of both raw and denoised spectra
print(result)
```

#### Denoising search on all spectra against reference library
```python
import spectral_denoising as sd
query_spectra= sd.read_msp('sample_data/query_spectra.msp')
reference_library =sd.read_msp('sample_data/reference_library.msp')

results = sd.denoising_search_batch(query_spectra['peaks'], query_spectra['precursor_mz'], reference_library) 
# results will be a list of all correspoinding precursor mz candidates, each one with entropy similarities of both raw and denoised spectra (using reference spectra melecular information)
print(results[0])# this will show denoising search result for the first spectra in msp file
```
## Working examples
More working examples can be found under notebooks directory.
