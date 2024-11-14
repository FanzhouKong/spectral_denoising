import spectral_denoising as sd
def main():
  query_spectra= sd.read_msp('sample_data/query_spectra.msp')
  reference_library =sd.read_msp('sample_data/reference_library.msp')
  
  results = sd.denoising_search_batch(query_spectra['peaks'], query_spectra['precursor_mz'], reference_library, )
  results = sd.denoising_search_batch(query_spectra['peaks'], query_spectra['precursor_mz'], reference_library, smiles_col = 'formula') # use this if you wish to provide formula information instead of smiles information
  # results will be a list of all correspoinding precursor mz candidates, each one with entropy similarities of both raw and denoised spectra (using reference spectra melecular information)
  print(results[0])# this will show denoising search result for the first spectra in msp file
if __name__ == "__main__":
    main()