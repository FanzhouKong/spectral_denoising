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

    Entropy similarity of spectrum with noise: 0.33.
    Entropy similarity of denoised spectrum: 1.00.


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

The very first step is to get the molecular formula of the precursor ion using the SMILES code. If the adduct contains atom other than a proton, the master formula will be modified accordingly to allow loss with the adduct. Moreover, 2 nitrogen and 1 oxygen is added if benzene substructure is present, to account for forming of rare adduct in the collision cells (`More info <https://chemistry-europe.onlinelibrary.wiley.com/doi/full/10.1002/cphc.201200313>`_.)

.. code-block:: python

    from spectral_denoising.spectral_denoising import *
    smiles = 'O=c1nc[nH]c2nc[nH]c12'
    adduct = '[M+Na]+'
    print(prep_formula(smiles, adduct))

The output will be:

.. code-block:: python

    'C5H4N4NaO'

Step 1: Get precursor ion infrmation
-----------------------------------------------------------
Then we want to get precursor statistics. If real precursor ion exist, the algorithm will prefer to use it since then the loss calculation will be free of systematic error. If not, the algorithm will use the computed precursor m/z. The real mass error would also be calculated for this step, and if it is larger than the mass_tolerance fed into the function, it will also be slightly increased to account for that.

For this reason, it is recommended to use a smaller mass tolerance to start with. 

At this step, the spectra is also sliced into 2 parts, the precursor region and fragment region, using ``slice_spectrum`` function. Since the algorithm focuses on the relative loss, only the fragment region is denoised, while the precursor region was kept intact and will be returned at the very last.

.. code-block:: python

    from spectral_denoising.spectral_denoising import *
    from spectral_denoising.chem_utils import *
    smiles = 'O=c1nc[nH]c2nc[nH]c12'
    adduct = '[M+Na]+'
    peak = np.array([[48.992496490478516 ,154.0],
                    [63.006099700927734, 265.0],
                    [79.02062225341797, 521.0],
                    [159.02373146795, 999]], 
                    
                    dtype = np.float32)
    computed_pmz = calculate_precursormz(adduct, smiles)
    pmz, mass_threshold = get_pmz_statistics(peak, computed_pmz, mass_tolerance=0.005)
    print(pmz, mass_threshold)

.. code-block:: python

    159.02373 0.006996155


Step2: Populate all possible subformulas from master formula
-------------------------------------------------------------

The next step is to populate all possible subformulas (with their masses) from the master formula. This can be easiily done with ``get_all_subfromulas`` function. The candidate formulas and masses are sorted so that search speed get facilitated.

.. code-block:: python

    all_possible_candidate_formula,all_possible_mass = get_all_subformulas(master_formula)
    print(all_possible_candidate_formula[1:5], all_possible_mass[1:5]) # only show first 5 for brevity

.. code-block:: python

    ['H', 'H2', 'H3', 'H4'] [1.00782503 2.01565006 3.0234751  4.03130013]


Step 3: Evaluate if a given ion could be formed from a plausible subformula loss
---------------------------------------------------------------------------------

For any given fragment ion, the algorithm will try to find a plausible subformula loss that could form this ion (function ``check_cnadidates``).

If such loss can be found, this ion will be given a tag 'True', otherwise a 'False' tag. This is done through function ``get_denoise_tag``.


Step 4: Retaining True fragment ions and add back the precursor ions
--------------------------------------------------------------------

Once the denoised tag was created, only ions with 'True' tag will be kept. The return spectra will be denoised fragment ions and precursor ions (function ``add_spectra``).

References
----------


.. autofunction:: spectral_denoising.spectral_denoising.formula_denoising
    :noindex:

.. autofunction:: spectral_denoising.spectral_denoising.prep_formula
    :noindex:

.. autofunction:: spectral_denoising.spectral_denoising.get_pmz_statistics
    :noindex:

.. autofunction:: spectral_denoising.spectral_denoising.get_all_subformulas
    :noindex:  

.. autofunction:: spectral_denoising.spectral_denoising.get_all_subformulas
    :noindex:

.. autofunction:: spectral_denoising.spectral_denoising.get_denoise_tag
    :noindex:

.. autofunction:: spectral_denoising.spectral_denoising.check_candidates
    :noindex:

.. autofunction:: spectral_denoising.seven_golden_rules.check_ratio
    :noindex:

.. autofunction:: spectral_denoising.spectral_operations.slice_spectrum
    :noindex:

.. autofunction:: spectral_denoising.spectral_operations.add_spectra
    :noindex:





