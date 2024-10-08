===========
Quickstart
===========

Once you have installed the package, you can quickly perform spectral denoising or denoising search on your data.
We also provided demo data `here <https://github.com/FanzhouKong/spectral_denoising/tree/main/sample_data>`_.

Spectral denoising
------------------

.. code-block:: python

    import spectral_denoising as sd
    query_data = sd.read_msp('sample_data/noisy_spectra.msp').iloc[0]
    query_spectrum, query_smiles, query_adduct, query_pmz = query_data['peaks'], query_data['smiles'], query_data['adduct'], query_data['precursor_mz']
    sd.spectral_denoising(query_spectrum, query_smiles, query_adduct, query_pmz)


Denoising search
----------------

.. code-block:: python

    import spectral_denoising as sd
    query_spectra= sd.read_msp('sample_data/query_spectra.msp')
    reference_library =sd.read_msp('sample_data/reference_library.msp')
    query_spectrum, query_pmz = query_spectra.iloc[0]['peaks'], query_spectra.iloc[0]['precursor_mz']
    sd.denoising_search(query_spectrum, query_pmz, reference_library)