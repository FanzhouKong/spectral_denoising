from toolsets.feature_alignment import find_feature
from tqdm import tqdm
import numpy as np
import pandas as pd
from toolsets.search import quick_search_values
from toolsets.flash_entropy_helper import flash_entropy
def identity_search_df(alignment, library,rt_offset = 5, use_rt =True):
    pl = alignment.dropna(subset=['msms'], ignore_index=True)
    pl['PeakID']=np.arange(len(pl))

    # pl = alignment.copy()
    library_sorted = library.sort_values(by = 'reference_precursor_mz', ascending=True)
    result = pd.DataFrame()
    pl_matched = pd.DataFrame()
    for index, row in tqdm(pl.iterrows(), total = len(pl)):
        if use_rt == True:
            candidate = find_feature(library, mz = row['pmz'], rt = row['rt'], mz_column='reference_precursor_mz', rt_column='rt', rt_offset=rt_offset)
        else:
            candidate = quick_search_values(library_sorted, value_start=row['pmz']-0.005, value_end=row['pmz']+0.005, ifsorted=True)
        if len(candidate)>0:

            entropy_temp = []
            for i in range(0, len(candidate)):
                simi = flash_entropy(peak1 = row['msms'], pmz1= row['pmz'], peak2 = candidate.iloc[i]['peaks'],pmz2=candidate.iloc[i]['reference_precursor_mz'])
                entropy_temp.append(simi['identity_search'][0])
            candidate.insert(0, 'peak_id', row['PeakID'])
            candidate.insert(1, 'identity_search_score', entropy_temp)
            candidate.insert(2, 'sample_pmz', row['pmz'])
            candidate.insert(3, 'sample_rt', row['rt'])
            candidate.insert(4, 'sample_peak', row['msms'])

            candidate = candidate[candidate['identity_search_score']>0.75]
            candidate.sort_values(by = 'identity_search_score', ascending=False, inplace=True)
            if len(candidate)>0:
                result = pd.concat([result, candidate], ignore_index=True, axis=0)
                row_df = pd.DataFrame([row]*len(candidate))
                # print(len(candidate))
                # return(row_df)
                row_df['annotation']=candidate['reference_name'].values
                row_df.insert(1, 'refrence_precursor_mz', candidate['reference_precursor_mz'].values)
                row_df.insert(3, 'reference_rt', candidate['rt'].values)
                row_df.insert(9, 'reference_msms', candidate['peaks'].values)
                row_df.insert(4, 'reference_smiles', candidate['reference_smiles'].values)
                row_df.insert(0, 'entropy_similarity', candidate['identity_search_score'].values)
                row_df.insert(6, 'BSD_ID', candidate['reference_mix'].values)
                pl_matched = pd.concat([pl_matched, row_df], ignore_index=True)
                # return(row_df, candidate)
    max_all = []
    for index, row in pl_matched.iterrows():
        max = row[row['reference_file']]
        max_all.append(max)
    pl_matched.insert(2, 'max_intensity', max_all)

    return (pl_matched)