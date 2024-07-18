import numpy as np
from tqdm import tqdm
from toolsets.search import quick_search_sorted
import toolsets.spectra_operations as so
def identity_search_single(msms, pmz, ref_lib, if_lib_sorted = True, remove_precursor = True, ms1_error = 0.01, ms2_error = 0.02):
    if if_lib_sorted == False:
        ref_lib.sort_values(by = 'PrecursorMZ', inplace = True, ascending=True)
    pmz_candidate = quick_search_sorted(ref_lib, 'PrecursorMZ', pmz -ms1_error, pmz+ms1_error)
    if len(pmz_candidate)==0:
        return(pmz_candidate)
    simi = []
    if remove_precursor == False:
        pmz = pmz+1.6
    for index, row in pmz_candidate.iterrows():
        simi.append(so.entropy_identity(msms, row['spectrum'], pmz, ms2_error = ms2_error))
    pmz_candidate['entropy_similarity']=simi
    pmz_candidate.sort_values(by = 'entropy_similarity', ascending = False, inplace = True)
    return pmz_candidate
    # row = pmz_candidate.iloc[0]
    # return([row['entropy_similarity'],row['spectrum'], row['SMILES'], row['Precursor_type'], row['Name']])
# def identity_search_all(df_input, ref_lib, if_lib_sorted = True, remove_precursor=True, ms1_error = 0.01, ms2_error = 0.02):
#     df = df_input.copy()
#     library_peaks = []
#     entropy_raw = []
#     reference_smiles =[]
#     reference_adduct = []
#     reference_name = []
#     if if_lib_sorted == False:
#         ref_lib.sort_values(by = 'PrecursorMZ', ascending=True, inplace = True)
#     for index, row in tqdm(df.iterrows(), total = len(df)):
#         temp_result=identity_search_single(row['peak'], row['precursor_mz'], ref_lib, if_lib_sorted = True, remove_precursor = remove_precursor,ms1_error = ms1_error,
#                                            ms2_error = ms2_error)
#         if len(temp_result)>0:
#             library_peaks.append(temp_result[1])
#             entropy_raw.append(temp_result[0])
#             reference_smiles.append(temp_result[2])
#             reference_adduct.append(temp_result[3])
#             reference_name.append(temp_result[4])
#         else:
#             library_peaks.append(np.NAN)
#             entropy_raw.append(np.NAN)
#             reference_smiles.append(np.NAN)
#             reference_adduct.append(np.NAN)
#             reference_name.append(np.NAN)
#     df['library_peaks']=library_peaks
#     df['entropy_raw']=entropy_raw
#     df['reference_name']=reference_name
#     df['reference_smiles']=reference_smiles
#     df['reference_adduct']=reference_adduct
#     return(df)