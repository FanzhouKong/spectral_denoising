import spectral_denoising as sd
query_spectra= sd.read_msp('sample_data/query_spectra.msp')
reference_library =sd.read_msp('sample_data/reference_library.msp')
query_spectrum, query_pmz = query_spectra.iloc[0]['peaks'], query_spectra.iloc[0]['precursor_mz'] # just the first spectrum
result = sd.denoising_search(query_spectrum, query_pmz, reference_library)
# result will return all precursor candidates of the query spectrum, each with entropy similarities of both raw and denoised spectra
print(result)