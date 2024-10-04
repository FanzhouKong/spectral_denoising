========================================
Spectral denoising: electronic denoising
========================================


The ``electronic_denoising`` function removes obvious electronic noise ions in MS/MS spectra, usually shown as multiple ions with identical intensities ( `Grass noise <https://mzmine.github.io/mzmine_documentation/module_docs/featdet_mass_detection/mass-detection.html>`_ )

According to empiracally tested on NIST23 database, in a given spectrum, the number of ions with identical intensities more than 4 is extremely unlikely (\< 0.05%). Thus, the ``electronic_denoising`` function removes ions with identical intensities greater than 4.

Here's an example:

.. code-block:: python

    import spectral_denoising as sd
    from spectral_denoising.noise import *
    peak = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
    pmz = 91
    noise = generate_noise(pmz, lamda=10, n = 50)
    peak_with_noise = add_noise(peak, noise)

    peak_denoised = sd.electronic_denoising(peak)
    print(f'Entropy similarity of spectra with noise: {sd.entropy_similairty(peak_with_noise,peak, pmz ):.2f}.')
    print(f'Entropy similarity of spectra with noise: {sd.entropy_similairty(peak_denoised,peak, pmz ):.2f}.')

The output will be:

.. code-block:: python

    Entropy similarity of spectra with noise: 0.37.
    Entropy similarity of spectra with noise: 1.00.

References
----------


.. autofunction:: spectral_denoising.electronic_denoising
    :noindex:

.. autofunction:: spectral_denoising.noise.generate_noise
    :noindex:

.. autofunction:: spectral_denoising.noise.generate_chemical_noise
    :noindex:

.. autofunction:: spectral_denoising.noise.add_noise
    :noindex:

