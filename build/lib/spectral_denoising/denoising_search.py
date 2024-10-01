from spectral_denoising.search_utils import quick_search_values
import pandas as pd
from tqdm import tqdm
import numpy as np
import multiprocessing as mp
from spectral_denoising.spectral_denoising import spectra_denoising
import spectral_denoising.spectral_operations as so
def denoising_search_batch(quene_msms, quene_pmz, reference_lib, identitiy_search_mass_error = 0.01, max_allowed_denoising_tolerance=0.005, 
                           pmz_col = 'precursor_mz', smiles_col = 'smiles',
                           adduct_col = 'adduct', msms_col = 'peaks'):
    if len(quene_msms)!= len(quene_pmz):
        print('The length of msms and pmz should be the same')
        return pd.DataFrame()

    with mp.Pool(processes=6) as pool:
            # Use starmap to handle multiple parameters
        results = pool.starmap(denoising_search, tqdm([(msms, pmz, reference_lib, identitiy_search_mass_error, max_allowed_denoising_tolerance, 
                            pmz_col, smiles_col, adduct_col,  msms_col) for msms, pmz in zip(quene_msms, quene_pmz)], total=len(quene_msms)))
    return results



def denoising_search(msms, pmz, reference_lib, identitiy_search_mass_error=0.01, max_allowed_denoising_tolerance = 0.005, 
                            pmz_col = 'precursor_mz', smiles_col = 'smiles', adduct_col = 'adduct', msms_col = 'peaks'):
    pmz_candidates = quick_search_values(reference_lib, pmz_col, pmz - identitiy_search_mass_error, pmz + identitiy_search_mass_error)
    
    if len(pmz_candidates) == 0:
        return pd.DataFrame()#no pmz matches
    for index, row in pmz_candidates.iterrows():
        try:
            msms_d = spectra_denoising(msms, row[smiles_col], row[adduct_col], max_allowed_deviation=max_allowed_denoising_tolerance)
        except:
            msms_d = np.nan
            continue
        pmz_candidates.loc[index, 'entropy_similarity'] = so.entropy_similairty(msms, row[msms_col], pmz = pmz, ms2_error = 0.01)
        pmz_candidates.loc[index, 'denoised_similarity'] = so.entropy_similairty(msms_d, row[msms_col], pmz = pmz, ms2_error = 0.01)
    pmz_candidates['quene_peaks']=[msms]*len(pmz_candidates)
    pmz_candidates['quene_pmz_denoised']=[msms_d]*len(pmz_candidates)
    pmz_candidates['quene_pmz']=pmz
    return pmz_candidates
