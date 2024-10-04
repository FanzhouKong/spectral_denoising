.. spectral_denoising documentation master file, created by
   sphinx-quickstart on Wed Oct  2 15:26:10 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

spectral_denoising documentation
================================

Welcome to python package for spectral denoising! The source code can be download from `Spectral denoising Github repository <https://github.com/FanzhouKong/spectral_denoising>`_.

If you encounter any issues, queries or need support, don't hesitate to contact :email:`Fanzhou Kong <fzkong@ucdavis.edu>`

For detail methodology or you want to cite this work, please refer to the paper: https://pubmed.ncbi.nlm.nih.gov/39108483/

This package contains the following modules:

1. **Classic spectral denoising:** These functions are used to remove noise ions in annotated spectra based on intensity distribution and subformula assignment for each ion.
2. **Denoising search:** In cases where the spectra are not annotated and molecular information are not feasible, spectral denoising is performed on query spectrum using molecular information of all candidate spectra within predefined precursor mass range. The denoised spectra are used to calculate similarity against reference spectra.

.. toctree::
   :maxdepth: 1
   :caption: Quickstart:

   install
   quickstart

.. toctree::
   :maxdepth: 1
   :caption: Spectral denoising
   
   sd_electronic_denoising
   sd_formula_denoising
   sd_single_spectrum_denoising
   sd_batch_spectra_denoising
   sd_useful_functions
   sd_api

.. toctree::
   :maxdepth: 1
   :caption: Denoising search
   
   ds_single_spectrum_search
   ds_batch_spectra_search
   ds_useful_functions
   ds_api
