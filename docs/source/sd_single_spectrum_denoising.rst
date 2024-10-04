=======================================
Spectral denoising: for single spectrum
=======================================

The ``spectral_denoising`` function, which is the core function of the project, removes noise ions in MS/MS spectra. 

The function essentially is a warpper function for performing ``electronic_denoising`` and ``formula_denoising`` functions in sequence. Just like ``formula_denoising``, ``spectral_denoising`` also requires molecular information (SMILES and adduct) to remove noise ions.
If no valid ion was left after electronic denoising, the function will return nan.

Example usage:

.. code-block:: python

    import spectral_denoising as sd
    from spectral_denoising.noise import *
    from spectral_denoising.spectral_operations import *
    from spectral_denoising.chem_utils import *
    peak = np.array([[48.992496490478516 ,154.0],
                  [63.006099700927734, 265.0],
                  [79.02062225341797, 521.0]], dtype = np.float32)
    noise = generate_noise(pmz, lamda=10, n = 50)
    smiles = 'O=c1nc[nH]c2nc[nH]c12'
    adduct = '[M+Na]+'
    pmz = calculate_precursormz(adduct,smiles)
    peak_with_noise = add_noise(peak, noise)
    peak_denoised = sd.spectral_denoising(peak_with_noise, smiles, adduct)
    print(f'Entropy similarity of spectra with noise: {sd.entropy_similairty(peak_with_noise,peak, pmz ):.2f}.')
    print(f'Entropy similarity of spectra with noise: {sd.entropy_similairty(peak_denoised,peak, pmz ):.2f}.')

The output will be:

.. code-block:: python

    Entropy similarity of spectra with noise: 0.47.
    Entropy similarity of spectra with noise: 0.99.

References
----------


.. autofunction:: spectral_denoising.spectral_denoising
    :noindex: