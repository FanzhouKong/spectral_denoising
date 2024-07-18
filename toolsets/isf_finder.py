import pandas as pd
from toolsets.search import string_search, num_search
import bisect
from tqdm import tqdm
from toolsets.spectra_operations import break_spectra, entropy_similarity_default
# def isf_finder_pandas(unknown, library):
#     isf_table = initialize_isf_table()
def isf_finder(msp_compound,msp_unknown, only_precursor_found = False):
    isf_table = initialize_isf_table()
    for index, row in tqdm(msp_unknown.iterrows(), total=len(msp_unknown)):
        candidates = num_search(msp_compound, 'RETENTIONTIME', float(row['RETENTIONTIME']),direction='between', step = 6/60, inclusion=True )
        entropy_similarity = []
        for candidacy_rt, candidate_rt in candidates.iterrows():

            msms_truncated = truncate_msms(msms = candidate_rt['peaks'], highest_allowed=row['PRECURSORMZ']+1.5)
            entropy_similarity.append(entropy_similarity_default(row['peaks'], msms_truncated))
        isf_parent_index = [i for i,v in enumerate(entropy_similarity) if v > 0.6]
        candidates = candidates.iloc[isf_parent_index]
        if len(candidates)>0:
            for candidacy, candidate in candidates.iterrows():
                    # print(candidacy)
                for column in ['NAME','PRECURSORMZ', 'PRECURSORTYPE','RETENTIONTIME','Comment','peaks']:
                    isf_table['isf'+'_'+column].append(row[column])
                for column in ['NAME','PRECURSORMZ', 'PRECURSORTYPE','RETENTIONTIME','Comment','peaks']:
                    isf_table['Compound'+'_'+column].append(candidate[column])
                msms_truncated = truncate_msms(msms = candidate['peaks'], highest_allowed=row['PRECURSORMZ']+1.5)
                isf_table['Entropy_similarity'].append(entropy_similarity_default(row['peaks'], msms_truncated))
                    # entropy_similarity.append(entropy_similarity_default(row['peaks'], msms_truncated))
                isf_table['Precursor_found'].append(find_precursor_in_ms2(candidate['peaks'], row['PRECURSORMZ']))    
            # break
            
    isf_table = pd.DataFrame.from_dict(isf_table)
    if only_precursor_found == True:
        isf_table = string_search(isf_table, 'Precursor_found', True)
    isf_table.reset_index(inplace = True, drop = True)
    return (isf_table)
def isf_finder_validation(msp_compound,msp_unknown, only_precursor_found = False):
    isf_table = initialize_isf_table()
    for index, row in tqdm(msp_unknown.iterrows(), total=len(msp_unknown)):
        candidates = num_search(msp_compound, 'RETENTIONTIME', float(row['RETENTIONTIME']),direction='between', step = 6/60, inclusion=True )
        if len(string_search(candidates, 'NAME', row['NAME']))>0:

            entropy_similarity = []

            for candidacy_rt, candidate_rt in candidates.iterrows():
                msms_truncated = truncate_msms(msms = candidate_rt['peaks'], highest_allowed=row['PRECURSORMZ']+1.5)
                entropy_similarity.append(entropy_similarity_default(row['peaks'], msms_truncated))
            isf_parent_index = [i for i,v in enumerate(entropy_similarity) if v > 0.6]
            candidates = candidates.iloc[isf_parent_index]
            if len(candidates)>0:
                for candidacy, candidate in candidates.iterrows():
                    # print(candidacy)
                    for column in ['NAME','PRECURSORMZ', 'PRECURSORTYPE','RETENTIONTIME','Comment','peaks']:
                        isf_table['isf'+'_'+column].append(row[column])
                    for column in ['NAME','PRECURSORMZ', 'PRECURSORTYPE','RETENTIONTIME','Comment','peaks']:
                        isf_table['Compound'+'_'+column].append(candidate[column])
                    msms_truncated = truncate_msms(msms = candidate['peaks'], highest_allowed=row['PRECURSORMZ']+1.5)
                    isf_table['Entropy_similarity'].append(entropy_similarity_default(row['peaks'], msms_truncated))
                    # entropy_similarity.append(entropy_similarity_default(row['peaks'], msms_truncated))
                    isf_table['Precursor_found'].append(find_precursor_in_ms2(candidate['peaks'], row['PRECURSORMZ']))
            # break
    isf_table = pd.DataFrame.from_dict(isf_table)
    if only_precursor_found == True:
        isf_table = string_search(isf_table, 'Precursor_found', True)
    isf_table.reset_index(inplace = True, drop = True)
    return (isf_table)
def initialize_isf_table():
    isf_table = {}
    for status in ['isf', 'Compound']:
        for column in ['NAME','PRECURSORMZ', 'PRECURSORTYPE','RETENTIONTIME','Comment','peaks']:
            isf_table[status+'_'+column] = []
    # isf_table['Compound_NAME']=[]
    isf_table['Entropy_similarity']=[]
    isf_table['Precursor_found'] =[]
    return(isf_table)
def find_precursor_in_ms2(msms, precursormz, step = 0.01):
    mass, intensity = break_spectra(msms)
    upper_allowed = bisect.bisect_right(mass, precursormz+step)
    lower_allowed = bisect.bisect_left(mass, precursormz-step)
    allowed_mass = mass[lower_allowed:upper_allowed]
    allowed_intensity = intensity[lower_allowed:upper_allowed]
    if len(allowed_intensity)>0 and len(allowed_mass)>0:
        if max(allowed_intensity)>0:
            return(True)
        else:
            return (False)
    else:
        return (False)

    # else
    # mass_precursor = mass[upper_allowed:len(mass)]
    # intensity_precursor = intensity[upper_allowed:len(intensity)]
    # if mass_precursor == []:
    #     mass_precursor.append(parention)
    #     intensity_precursor.append(-1)
