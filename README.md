[![Test Spectral denoising package](https://pypi.org/project/spectral-denoising/)]
If you have any questions, feel free to send me E-mails: fzkong@ucdavis.edu. If you find this package useful, please consider citing the following papers:

> . **Denoising Search doubles the number of metabolite and exposome annotations in human plasma using an Orbitrap Astral mass spectrometer**, Res Sq, 2024). [https://doi.org/10.1038/s41592-023-02012-9]


# Project information

Remove noise ions from MS/MS spectra has been tackled for years by mass spectrometrists. Noise ions in MS/MS spectra are largely categorized as 1. electronic noises and 2. chemical noises. 
In this project, we aim to eliminate both chemical noise and electronic noises for improving high-confidence compound identification. 
Integrating such process into spectra matching process, we developed denoising search, which psudo-denoise spectra based on molecular information fetched from reference databases.
This project also provides useful tools to read, write, visualize and compare spectra.

# How to use this package

This repository in Python.

### Installation

```bash
pip install spectral-denoising
```

### Usage of Classical spectral denoising (electronic denoising and chemical denoising)

```python
import numpy as np
import spectral_denoising as sd
import spectral_denoising.spectral_operations as so

quene_data = sd.read_msp('../sample_data/noisy_spectra.msp').iloc[0] # just use the first spectra
reference_data= sd.read_msp('../sample_data/clean_spectra.msp').iloc[0] 
quene_spectra, quene_smiles, quene_adduct, quene_pmz = quene_data['peaks'], quene_data['smiles'], quene_data['adduct'], quene_data['precursor_mz']
reference_spectra = reference_data['peaks']

# Visualize the head-to-tail plots
sd.head_to_tail_plot(quene_spectra, reference_spectra, pmz = quene_pmz)

# Run spectra denoising will remove electronic noises and chemical noises.
msms_denoised = sd.spectra_denoising(quene_spectra, quene_smiles, quene_adduct) #denoise the spectrum based on the smiles/adduct information

# Calculate spectral entropy before and after denoising.
print(f'the raw spectra has spectral entropy of {so.spctrum_entropy(quene_spectra):.2f}')
print(f'the denoised spectra has spectral entropy of {so.spctrum_entropy(msms_denoised):.2f}')
# Calculate entropy similarity.
sd.head_to_tail_plot(msms_denoised, reference_spectra, pmz = quene_pmz)
```
### Example of usage for batch mode can be found in /notebook/spectral_denoising_demo.ipynb


### Usage of Flash Entropy Search

```python
import spectral_denoising as sd
data_dir = '../sample_data/'
quene_spectra= sd.read_msp(os.path.join(data_dir, 'quene_spectra.msp'))
reference_library =sd.read_msp(os.path.join(data_dir, 'sample_library.msp'))
quene_spectrum, quene_pmz = quene_spectra.iloc[0]['peaks'], quene_spectra.iloc[0]['precursor_mz']
sd.denoising_search(quene_spectrum, quene_pmz, reference_library)
```

### Example of usage for batch mode can be found in /notebook/denoising_search_demo.ipynb


Also, refer to [MSViewer repository](https://github.com/YuanyueLi/MSViewer) for a working example of using this package in a web application.
