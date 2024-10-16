import pandas as pd
import numpy as np
import sys
import inspect
import os
sys.path.insert(0, '..') 
import spectral_denoising as sd
def main():
	query_spectra= sd.read_df('/Users/fanzhoukong/Documents/GitHub/Libgen_data/Astral/Astral/alignment/pos/astral.csv')

	reference_lib= sd.read_df('/Users/fanzhoukong/Documents/GitHub/Libgen_data/curated_library/csv/pos_orbi_sorted.csv')
	reference_lib.dropna(subset=['smiles'], inplace=True)
	query_spectra = query_spectra[['precursor_mz','peaks']]
	query_spectra.dropna(subset=['peaks'],inplace=True)
	results = sd.denoising_search_batch(query_spectra['peaks'], query_spectra['precursor_mz'], reference_lib)
if __name__=="__main__":
	main()
