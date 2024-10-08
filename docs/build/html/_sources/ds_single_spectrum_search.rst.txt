========================================
Denoising search: for single spectrum
========================================

The ``denoising_search`` function, which is the core function of the project, performs identity search with spectral denoising integrated.
In this way, the query spectrum is denoised using all molecular information obtained from candidate spectra with predefined precursor mass range. The entropy similarity scores are computed with denoised spectra.

Example usage:
The demo data can be found `here <https://github.com/FanzhouKong/spectral_denoising/tree/main/sample_data>`_.

.. code-block:: python

    import pandas as pd
    import spectral_denoising as sd
    quene_spectra= sd.read_msp('sample_data/query_spectra.msp')
    reference_library =sd.read_msp('sample_data/reference_library.msp')
    quene_spectrum, quene_pmz = quene_spectra.iloc[0]['peaks'], quene_spectra.iloc[0]['precursor_mz']
    result = sd.denoising_search(quene_spectrum, quene_pmz, reference_library)

The result will give all candidate spectra within the precursor mass range, with additional column of 'query_peaks' (query spectrum), 'query_peaks_denoised' (denoised query spectra), 'entrpy_similarity' (entropy similarity of query spectra to reference spectra), and 'denoised_similarity' (entropy similarity of denoised query spectra to reference spectra).

References
----------

.. autofunction:: spectral_denoising.denoising_search
    :noindex: