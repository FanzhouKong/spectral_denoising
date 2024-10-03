# Spectral denoising and denoising search

Test the package at: https://pypi.org/project/spectral-denoising/


If you have any questions, feel free to send me E-mails: fzkong@ucdavis.edu. If you find this package useful, please consider citing the following papers:

> . **Denoising Search doubles the number of metabolite and exposome annotations in human plasma using an Orbitrap Astral mass spectrometer**, Res Sq, 2024). [https://doi.org/10.1038/s41592-023-02012-9]


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

## Simple usage on single spectra
```python
import numpy as np
import spectral_denoising as sd
from spectral_denoising.spectral_operations import *
from spectral_denoising.chem_utils import *

smiles = 'O=c1nc[nH]c2nc[nH]c12'
adduct = '[M+Na]+'
pmz = calculate_precursormz(adduct,smiles)
peak = np.array([[48.992496490478516 ,154.0],
                  [63.006099700927734, 265.0],
                  [79.02062225341797, 521.0]], dtype = np.float32)
print(f'the spectrum entropy is {spctrum_entropy(peak):.2f}, the normalized entropy is {normalized_entropy(peak):.2f}')
# alternatively, you can store mass and intensity in separate arrays, and use pack_spectrum(mass, intensity) to get the peaks array
# e.g.mass,intensity = [48.992496490478516, 63.006099700927734, 79.02062225341797], [154.0, 265.0, 521.0]
# peak = pack_spectrum(mass, intensity)

# generate some electronic noise and add it to the peaks
from spectral_denoising.noise import *
noise = generate_noise(pmz, lamda=10, n = 50)
peak_with_noise = add_noise(peak, noise)
# use head_to_tail_plot to visualize the spectra, only in jupyter notebook
# sd.head_to_tail_plot(peaks_with_noise,peaks ,pmz)
print(f'the spectrum entropy is {spctrum_entropy(peak_with_noise):.2f}, the normalized entropy is {normalized_entropy(peak_with_noise):.2f}')
print(f'the entropy similarity of contaminated spectrum and the raw spectrum is {entropy_similairty(peak_with_noise,peak,  pmz = pmz):.2f}')

# perform spectral denosing and compare against the raw spectrum
peak_denoised = sd.spectral_denoising(peak_with_noise, smiles, adduct)
print(f'the entropy similarity of denoised spectrum and the raw spectrum is {entropy_similairty(peak_denoised, peak, pmz = pmz):.2f}')
# use head_to_tail_plot to visualize the spectra, only in jupyter notebook
# sd.head_to_tail_plot(peaks_denoised,peaks ,pmz)
```
### Spectral denoising on the all spectra from whole .msp file
```python
import spectral_denoising as sd
query_data = sd.read_msp('sample_data/noisy_spectra.msp')
query_peaks,query_smiles,query_adduct, query_pmz = query_data['peaks'],query_data['smiles'],query_data['adduct'], query_data['precursor_mz']
desnoied_peaks = sd.spectra_denoising_batch(quene_peaks,quene_smiles,quene_adduct) # this will return all denoised spectra in a list
```

### Usage of Denoising search

```python
import spectral_denoising as sd
data_dir = '../sample_data/'
quene_spectra= sd.read_msp(os.path.join(data_dir, 'quene_spectra.msp'))
reference_library =sd.read_msp(os.path.join(data_dir, 'sample_library.msp'))
quene_spectrum, quene_pmz = quene_spectra.iloc[0]['peaks'], quene_spectra.iloc[0]['precursor_mz']
sd.denoising_search(quene_spectrum, quene_pmz, reference_library)
```

### More example of usage for batch mode can be found in /notebook/denoising_search_demo.ipynb
