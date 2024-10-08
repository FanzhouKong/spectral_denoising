==============================
Denoising search: bacth mode
==============================

The ``denoising_search_batch`` function is essentially a wrapper function of ``denoising_search`` for batch data, while implemented in parallel processing. 
The function takes similari parameters as ``denoising_search``, but msms, pmz are now list or iteratable objects instead of single spectrum and float.

Example usage:
The demo data can be found `here <https://github.com/FanzhouKong/spectral_denoising/tree/main/sample_data>`_.

.. code-block:: python

    import pandas as pd
    import spectral_denoising as sd
    quene_spectra= sd.read_msp('sample_data/query_spectra.msp')
    reference_library =sd.read_msp('sample_data/reference_library.msp')
    results = sd.denoising_search_batch(quene_spectra['peaks'], quene_spectra['precursor_mz'],reference_library)

The results will be a list. At each index, it will give all candidate spectra with denoised information, just as in ``denoising_search``.

References
----------


.. autofunction:: spectral_denoising.denoising_search_batch
    :noindex: