=====================================
Spectral denoising: formula denoising
=====================================

The ``formula_denoising`` function removes chemical noise ions in MS/MS spectra by evaluating if it could be formed from a chemically plausible subformula loss to the precursor ion.

Basic usage
============

The ``formula_denoising`` function is used to remove noise ions in an annotated spectra. Thus, the molecular information is needed for this function (SMILES and adduct).

The mass tolerance is used to search the subformula loss of a given fragment ion. It is recommended to start with a smaller value. If the actual mass difference is larger than the mass tolerance (determined by precursor ion), this value will be increased. The max mass range to search for precursor ion is +/- 10 mDa.

The formula_db can be found at: `Formula_db <https://drive.google.com/file/d/1pEXiGc5l0YjRGfCEZXW7-Wz6D1dOSBxA/view?usp=drive_link>`_. It is already sorted by mass.

.. code-block:: python
    
        import spectral_denoising as sd
        from spectral_denoising.noise import *
        from spectral_denoising.chem_utils import *
        peak = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        smiles = 'C1=CC=CC=C1'
        adduct = '[M+H]+'
        pmz = calculate_precursormz(smiles, adduct)
        noise = generate_chemical_noise(pmz, lamda=10,polarity='+', formula_db, n = 10)
        peak_with_noise = add_noise(peak, noise)
    
        peak_denoised = sd.formula_denoising(peak, 'C1=CC=CC=C1', '[M+H]+')
        print(f'Entropy similarity of spectra with noise: {sd.entropy_similairty(peak_with_noise,peak, pmz ):.2f}.')
        print(f'Entropy similarity of denoised spectra: {sd.entropy_similairty(peak_denoised,peak, pmz ):.2f}.')

The output will be:

.. code-block:: python

    Entropy similarity of spectra with noise: 0.33.
    Entropy similarity of spectra with noise: 1.00.


References
----------


.. autofunction:: spectral_denoising.formula_denoising
    :noindex:

.. autofunction:: spectral_denoising.noise.generate_chemical_noise
    :noindex:
    

Want to know details about implementation?
===========================================



Step 0: Modify the master formula based on SMILES and adduct information
-------------------------------------------------------------------------

Step 1: Get precursor ion infrmation
-----------------------------------------------------------

Step2: Populate all possible subformulas from master formula
-------------------------------------------------------------

Step 3: Evaluate if a given ion could be formed from a plausible subformula loss
---------------------------------------------------------------------------------

Step 4: Filter out ions that could not be formed from a plausible subformula loss
----------------------------------------------------------------------------------

Step 5: Add back the precursor ions
------------------------------------




