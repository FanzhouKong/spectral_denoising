=================================
Spectral denoising: in bacth mode
=================================


The ``spectral_denoising_batch`` is just a simple wrapper function that performs ``spectral_denoising`` in batch mode, with parallelization implementaion. 

The sample data can be found `here <https://github.com/FanzhouKong/spectral_denoising/tree/main/sample_data>`_.

Example usage:

.. code-block:: python

    import spectral_denoising as sd
    from spectral_denoising.file_io import *
    query_data = sd.read_msp('sample_data/noisy_spectra.msp')
    query_spectra,query_smiles,query_adducts = query_data['peaks'],query_data['smiles'],query_data['adduct']
    denoised_spectra = sd.spectral_denoising_batch(query_spectra,query_smiles,query_adducts)


The output 'denoised_spectra' will be a list of denoised spectra, with the same order as the input spectra. 

References
----------


.. autofunction:: spectral_denoising.spectral_denoising_batch
    :noindex:

.. autofunction:: spectral_denoising.file_io.read_msp
    :noindex: