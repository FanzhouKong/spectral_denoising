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

This repository in Python. A python version >= 3.8 is preferred, and must be < 3.13.

Detailed documentation can be found at: https://spectral-denoising.readthedocs.io/en/latest/index.html

### Installation

```bash
pip install spectral-denoising
```

### Usage of Classic spectral denoising (electronic denoising and chemical denoising)
The demo data used here can be found under sample_data directory.

**Note: Even all functions have a default 'smiles' information column, the function would also accept formula as input. If wanted, just replace the the smiles with formula information (for single run) or replace 
column name of 'smiles' to your column name that contains formula information (batch mode).**

#### Simple usage on single spectra
**Note: if you try to use the batch mode in script and compile it in terminal, please wrap the code in main() function since they are implemented in parallal with multiprocessing and directly calling it will cause issues.**
```python
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
```
#### Spectral denoising on the all spectra from .msp file
```python
import spectral_denoising as sd
def main():
  query_data = sd.read_msp('sample_data/noisy_spectra.msp')
  query_peaks,query_smiles,query_formula,query_adduct, query_pmz = query_data['peaks'],query_data['smiles'],query_data['formula'],query_data['adduct'], query_data['precursor_mz']
  desnoied_peaks = sd.spectral_denoising_batch(query_peaks,query_smiles,query_adduct) # this will return all denoised spectra in a list
  # desnoied_peaks = sd.spectral_denoising_batch(query_peaks,query_formula,query_adduct) # this two would both work, if you wish to use formula information instead of SMILES information
if __name__ == "__main__":
    main()
```

### Usage of Denoising search
#### Denoising search on a single spectrum against reference library
```python
import spectral_denoising as sd
query_spectra= sd.read_msp('sample_data/query_spectra.msp')
reference_library =sd.read_msp('sample_data/reference_library.msp')
query_spectrum, query_pmz = query_spectra.iloc[0]['peaks'], query_spectra.iloc[0]['precursor_mz'] # just the first spectrum
result = sd.denoising_search(query_spectrum, query_pmz, reference_library) # default it will use SMILES information from 'smiles' column
# result = sd.denoising_search(query_spectrum, query_pmz, reference_library, smiles_col = 'formula') # use this if you wish to provide formula information instead of smiles information
# result will return all precursor candidates of the query spectrum, each with entropy similarities of both raw and denoised spectra
print(result)
```

#### Denoising search on all spectra against reference library
```python
import spectral_denoising as sd
def main():
  query_spectra= sd.read_msp('sample_data/query_spectra.msp')
  reference_library =sd.read_msp('sample_data/reference_library.msp')
  
  results = sd.denoising_search_batch(query_spectra['peaks'], query_spectra['precursor_mz'], reference_library, )
  # results = sd.denoising_search_batch(query_spectra['peaks'], query_spectra['precursor_mz'], reference_library, smiles_col = 'formula') # use this if you wish to provide formula information instead of smiles information
  # results will be a list of all correspoinding precursor mz candidates, each one with entropy similarities of both raw and denoised spectra (using reference spectra melecular information)
  print(results[0])# this will show denoising search result for the first spectra in msp file
if __name__ == "__main__":
    main()
```
Code for quick starts can be found in script director, with 2 demo files: spectral_denoising_demo.py and denoising_search_demo.py. Directly compiling these 2 files will produce similar results
## Working examples and reproducing results for the publication
All necessory data can be found here: https://drive.google.com/drive/folders/1xSKtLqNXukj6V8qP9c_e5MAbmbDLg2Ml?dmr=1&ec=wgc-drive-hero-goto

