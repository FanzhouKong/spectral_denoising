from .search_utils import quick_search_values
import pandas as pd
from tqdm import tqdm
import numpy as np
import multiprocessing as mp
from .spectral_denoising import spectral_denoising
from .spectral_operations import entropy_similairty
from .file_io import standardize_col
def denoising_search_batch(msms_query, pmz_query, reference_lib, identitiy_search_mass_error = 0.01, mass_tolernace=0.005,
                           pmz_col = 'precursor_mz', 
                           smiles_col = 'smiles',
                           adduct_col = 'adduct', msms_col = 'peaks'):
    """
    Perform batch denoising search on given MS/MS data and precursor m/z values with parallel processing.

    Parameters:
        msms_query (list): List of MS/MS spectra to be denoised.

        pmz_query (list): List of precursor m/z values corresponding to the MS/MS spectra.

        reference_lib (pandas.DataFrame): Reference library containing known spectra for comparison.

        identitiy_search_mass_error (float, optional): Mass error tolerance for identity search. Default is 0.01.

        mass_tolerance (float, optional): Maximum allowed tolerance for denoising. Default is 0.005.

        pmz_col (str, optional): Column name for precursor m/z in the reference library. Default is 'precursor_mz'.

        smiles_col (str, optional): Column name for SMILES in the reference library. Default is 'smiles'.

        adduct_col (str, optional): Column name for adducts in the reference library. Default is 'adduct'.

        msms_col (str, optional): Column name for MS/MS peaks in the reference library. Default is 'peaks'.
    Returns:
        pandas.DataFrame: DataFrame containing the results of the denoising search. Each index in the result DataFrame corresponds to the denoising search result of the corresponding input MS/MS spectrum.
    """
    print('ttt')
    if len(msms_query)!= len(pmz_query):
        print('The length of msms and pmz should be the same')
        return pd.DataFrame()

    with mp.Pool(processes=6) as pool:
            # Use starmap to handle multiple parameters
        results = pool.starmap(denoising_search, tqdm([(msms, pmz, reference_lib, identitiy_search_mass_error, mass_tolernace, 
                            pmz_col, smiles_col, adduct_col,  msms_col) for msms, pmz in zip(msms_query, pmz_query)], total=len(pmz_query)))
    
    return results



def denoising_search(msms, pmz, reference_lib, identitiy_search_mass_error=0.01, mass_tolernace = 0.005, 
                     pmz_col = 'precursor_mz', smiles_col = 'smiles', adduct_col = 'adduct', msms_col = 'peaks'):
    """
    Perform a denoising search on mass spectrometry data.
    In addition to traditional identity search, the query msms will be denoised using molecular information of all candidates within predefined identitiy_search_mass_error range.
    Then, entropy similarity scores are calculated for both the original and denoised spectra.

    Parameters:
        msms (np.array): The mass spectrometry data (peaks) to be denoised.

        pmz (float): The precursor m/z value associated with the query MS/MS spectra.

        reference_lib (pandas.DataFrame): The reference library containing known spectra and associated metadata.

        identitiy_search_mass_error (float, optional): The mass error tolerance for the identity search. Default is 0.01.

        mass_tolernace (float, optional): The mass tolerance for spectral denoising. Default is 0.005.

        pmz_col (str, optional): The column name for precursor m/z values in the reference library. Default is 'precursor_mz'.

        smiles_col (str, optional): The column name for SMILES strings in the reference library. Default is 'smiles'.

        adduct_col (str, optional): The column name for adduct information in the reference library. Default is 'adduct'.

        msms_col (str, optional): The column name for MS/MS peaks in the reference library. Default is 'peaks'.
    Returns:
        pd.DataFrame: A DataFrame containing all the candidate spectra, with entropy similarity scores calculated based on both denoised and raw spectra.
    """
    pmz_candidates = quick_search_values(reference_lib, pmz_col, pmz - identitiy_search_mass_error, pmz + identitiy_search_mass_error)
    
    if len(pmz_candidates) == 0:
        return pd.DataFrame()#no pmz matches
    for index, row in pmz_candidates.iterrows():
        try:
            msms_d = spectral_denoising(msms, row[smiles_col], row[adduct_col], mass_tolerance=mass_tolernace)
        except:
            msms_d = np.nan
            continue
        pmz_candidates.loc[index, 'entropy_similarity'] = entropy_similairty(msms, row[msms_col], pmz = pmz, ms2_error = 0.01)
        pmz_candidates.loc[index, 'denoised_similarity'] = entropy_similairty(msms_d, row[msms_col], pmz = pmz, ms2_error = 0.01)
    pmz_candidates['query_peaks']=[msms]*len(pmz_candidates)
    pmz_candidates['query_peaks_denoised']=[msms_d]*len(pmz_candidates)
    pmz_candidates['query_pmz']=pmz
    cols = ['name', 'cas','adduct','smiles','query_pmz', 'query_peaks','query_peaks_denoised','entropy_similarity','denoised_similarity']
    df_return = pmz_candidates[cols]

    df_return.sort_values(by='denoised_similarity', ascending=False, inplace=True)
    df_return = standardize_col(df_return)
    return df_return
