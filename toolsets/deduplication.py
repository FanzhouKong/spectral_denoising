import pandas as pd
# import toolsets.conversion_function as cf
from tqdm import tqdm
from molmass import Formula
import numpy as np
from toolsets.spectra_operations import weighted_average_spectra
from toolsets.search import string_search, num_search, quick_search_sorted
import toolsets.spectra_operations as so
def deduplication(data, reference_col = 'reference_inchikey', rt_offset = 5/60):
    matched_refined = pd.DataFrame()
    for key in tqdm(data[reference_col].unique()):
        data_temp = string_search(data, reference_col, key)
        master_seed = data_temp.iloc[np.argmax(data_temp['peak_apex_intensity'])]
        data_retain = quick_search_sorted(data_temp, 'rt', master_seed['rt'], step=rt_offset, ifsorted=False)
        data_retained=add_comment_retained(data_retain)
        data_recheck = data_temp[~data_temp['library_id'].isin(data_retained['library_id'])]
        data_recheked=recheck(data_retained, data_recheck)
        matched_refined = pd.concat([matched_refined, data_retained, data_recheked], ignore_index=True)
    return(matched_refined)
        # seed_col =
    # pass
def recheck(data_retained, data_recheck):
    comments = []
    updated_names = []
    rechecked = pd.DataFrame()
    for adduct in data_recheck['reference_adduct'].unique():
        data_retained_adduct = string_search(data_retained, 'reference_adduct', adduct)
        data_recheck_adduct = string_search(data_recheck, 'reference_adduct', adduct)
        if len(data_retained_adduct)>0:
            threshold = get_seed_data(data_retained_adduct,'peak_apex_intensity')/3
            reference_peak = get_seed_data(data_retained_adduct,'peaks')
            # print(threshold)
            for index, row in data_recheck_adduct.iterrows():
                entropy_similarity = so.entropy_similarity_default(reference_peak, row['peaks'])
                # print(entropy_similarity)
                if row['peak_apex_intensity']>=threshold and entropy_similarity >0.7:
                    comments.append(str(np.round(entropy_similarity,2))+'/'+str(np.round(row['peak_apex_intensity']/(threshold*3)*100, 1)))
                    updated_names.append(row['reference_name']+'_minor')
                    # print(threshold*3)
                    rechecked = rechecked.append(row)
    rechecked['comments']=comments
    rechecked['reference_name']=updated_names
    return(rechecked)
def add_comment_retained(data_retained):
    data_retained_return = pd.DataFrame()
    updated_name = []
    comments = []
    for adduct in data_retained['reference_adduct'].unique():
        data_temp = string_search(data_retained, 'reference_adduct', adduct)
        master_seed = data_temp.iloc[np.argmax(data_temp['peak_apex_intensity'])]
        for index, row in data_temp.iterrows():
            entropy_similarity=so.entropy_similarity_default(master_seed['peaks'],
                                                             row['peaks']
                                                             )
            comments.append(str(np.round(entropy_similarity,2))+'/'+str(np.round(row['peak_apex_intensity']/master_seed['peak_apex_intensity']*100,1)))
            if row['peak_apex_intensity']==master_seed['peak_apex_intensity']:
                updated_name.append(row['reference_name']+'_major')
            else:
                updated_name.append(row['reference_name']+'_'+str(np.round(row['peak_apex_intensity']/master_seed['peak_apex_intensity']*100,1))+"%")
            # print(row['reference_name'])
        data_retained_return=pd.concat([data_retained_return, data_temp], ignore_index=True)
    data_retained_return['reference_name']=updated_name
    data_retained_return['comments']=comments
    return(data_retained_return)
def get_seed_data(data_temp, return_col,seed_col = 'peak_apex_intensity'):
    master_seed = data_temp.iloc[np.argmax(data_temp[seed_col])]
    return(master_seed[return_col])

#
#
def match_by_rt(data_all, rt_tolerance=5/60, matching_column = 'reference_inchikey'):
    print("finding seed...")
    for inchi in tqdm(data_all[matching_column].unique()):
        data_temp = string_search(data_all, matching_column,inchi14)
        data_temp['seed_found']=_find_seed_rt(data_temp, rt_tolerance=rt_tolerance, iterative=True)
        matched_rt = matched_rt.append(data_temp)
    matched_rt_found = string_search(matched_rt, 'seed_found',-1, reverse=True)
    matched_rt_not_found = string_search(matched_rt, 'seed_found',-1)
    matched_rt_mapped = pd.DataFrame()
    print("making matches")
    for inchi14 in tqdm(matched_rt_found[matching_column].unique()):
        data_temp = string_search(matched_rt_found, matching_column, inchi14)
        if data_temp.iloc[0]['seed_found']!=-1:
            rt_seed = data_temp.iloc[0]['seed_found']
            rt_temp = data_temp.loc[rt_seed]['retention_time']
            data_temp_mapped = num_search(data_temp, 'retention_time',number = rt_temp, direction="between", step=10/60, inclusion=True)
            data_temp_mapped.drop_duplicates('rt_match_id', keep="first", inplace=True)
            matched_rt_mapped =  matched_rt_mapped.append(data_temp_mapped)
    matched_rt_unmapped = matched_rt_found[~matched_rt_found['rt_match_id'].isin(matched_rt_mapped['rt_match_id'])]
    for data_temp in [matched_rt_mapped, matched_rt_unmapped, matched_rt_not_found]:
        data_temp.reset_index(inplace= True, drop = True)
    matched_rt = []
    for inchi in matched_rt_mapped[matching_column].unique():
        compound_temp = string_search(matched_rt_mapped, matching_column, inchi)
        matched_rt.extend([string_search(compound_temp, 'rt_match_id',compound_temp.iloc[0]['seed_found']).iloc[0]['retention_time']]*len(compound_temp))
    matched_rt_mapped['matched_rt']=matched_rt
    return(matched_rt_mapped, matched_rt_unmapped, matched_rt_not_found)

def _find_seed_rt(data_temp, rt_hits_column = 'hits_rt', rt_tolerance = 5/60, iterative = False):

    if rt_hits_column not in data_temp:
        data_temp = _find_hit_rt(data_temp, rt_tolerance=rt_tolerance)
    if len(data_temp)<=1:
        return(data_temp.iloc[0]['library_id'])
        # print("the rt hits have not been found yet")
        # return(None)
    rt_hit_max = data_temp[rt_hits_column].max()
    if rt_hit_max ==0 and len(data_temp)>1:
        if iterative ==False:
            return(-1)
        else:
            data_temp.sort_values(by = 'peak_apex_intensity', inplace=True, ascending = False)
            return(data_temp.iloc[0]['library_id'])
    elif rt_hit_max ==0 and len(data_temp)==1:
        return(data_temp['peak_apex_intensity'].idxmax())
    seed_candidate = num_search(data_temp,rt_hits_column, rt_hit_max, direction='==')
    # return(seed_candidate)
    return(seed_candidate.loc[seed_candidate['peak_apex_intensity'].idxmax()]['library_id'])

#
#
#
def _find_hit_rt(data_temp, retention_time = 'rt', rt_tolerance=5/60):
    return_df = pd.DataFrame()
    hits_rt = []
    for index, row in data_temp.iterrows():
        rt_match_temp = num_search(data_temp, retention_time,number = row[retention_time], direction="between", step=rt_tolerance, inclusion=True)
        hits_rt.append(len(rt_match_temp)-1)
    data_temp['hits_rt']=hits_rt
    return(data_temp)
#
#     # for adduct in adducts:
#     #     hits = []
#     #     data_adduct_temp = string_search(data_temp, 'reference_adduct', adduct)
#     #     # data_adduct_temp_r = string_search(data_temp, 'reference_adduct', adduct, reverse=True)
#     #     for index, row in data_adduct_temp.iterrows():
#     #         rt_match_temp = num_search(data_adduct_temp_r, retention_time,number = row[retention_time], direction="between", step=6/60, inclusion=True)
#     #         hits.append(len(rt_match_temp))
#     #     data_adduct_temp['hits_rt']=hits
#     #     return_df = return_df.append(data_adduct_temp)



# def split_unique_duplicates(data, key_column = "key"):
#     print("i am in new")
#     # data = num_search(data, "ms1_intensity_ratio", 0.7, direction = ">", inclusion = True)
#     matched_confirmed = data.drop_duplicates('key', keep=False)
#     temp_idx = data.index.isin(matched_confirmed.index)
#     matched_duplicates = data[~temp_idx]
#     matched_confirmed.reset_index(inplace = True, drop = True)
#     # matched_confirmed=num_search(matched_confirmed, "ms1_intensity_ratio", 0.1, direction = ">", inclusion = True)
#     # peaks_cleaned = []
#     # for i in range(len(matched_confirmed)):
#     #     peaks_cleaned.append(weighted_average_spectra(matched_confirmed[i:i+1]))
#     # matched_confirmed['peaks']=peaks_cleaned
#     matched_duplicates.reset_index(inplace = True, drop = True)
#     # matched_duplicates=num_search(matched_duplicates, "ms1_intensity_ratio", 0.1, direction = ">", inclusion = True)
#     return(matched_confirmed, matched_duplicates)
# def match_by_rt(matched_confirmed, matched_duplicates, rt_tolerance = 6/60, referencing_column = 'reference_inchikey', rt_column = 'retention_time'):
#     common_compounds = list(set(matched_confirmed.reference_inchikey) & set(matched_duplicates.reference_inchikey))
#     # common_compounds = list(set(matched_confirmed['reference_Substrate Name']) & set(matched_duplicates['reference_Substrate Name']))
#     rt_matched = pd.DataFrame()
#     # i = 0
#     for inchi in tqdm(common_compounds):
#         # print(" i am matching ", i)
#         # i = i+1
#         unmatched_temp = string_search(matched_duplicates, referencing_column, inchi)
#         matched_temp = string_search(matched_confirmed, referencing_column, inchi)
#         for index, row in matched_temp.iterrows():
#             rt_confirmed = row[rt_column]
#         # rt_confirmed = string_search(matched_confirmed, "reference_InChIKey", inchi).iloc[0]['retention_time']
#             unmatched_candidates = num_search(unmatched_temp, rt_column, rt_confirmed, direction="between",
#                                               step = rt_tolerance, inclusion=True)
#             if len(unmatched_candidates)!=0:
#                 for adduct in unmatched_candidates['reference_adduct'].unique():
#                     unmatched_candidates_adduct = string_search(unmatched_candidates, "reference_adduct", adduct)
#                     # wa_spec = weighted_average_spectra(unmatched_candidates_adduct)
#                     # unmatched_candidates_adduct.at[0, 'peaks']=weighted_average_spectra(unmatched_candidates_adduct)
#                     # unmatched_candidates_adduct.iloc[0, unmatched_candidates_adduct.columns.get_loc('peaks')] = weighted_average_spectra(unmatched_candidates_adduct)
#                     rt_matched = pd.concat([rt_matched, unmatched_candidates_adduct], ignore_index=True)
#                 break
            
#     matched_confirmed_addi_matched = pd.concat([matched_confirmed, rt_matched], axis = 0, ignore_index = True)
#     matched_leftover = matched_duplicates[~matched_duplicates[referencing_column].isin(matched_confirmed_addi_matched[referencing_column])]
#     return(matched_confirmed_addi_matched, matched_leftover)
def check_for_missing_compounds( std_list,matched, check_column_std, check_column_matched):
    missing_compound = list(set(std_list[check_column_std]) - set(matched[check_column_matched]))
    missing_compounds_list = std_list.loc[std_list[check_column_std].isin(missing_compound)]
    return(missing_compounds_list)


# def deduplicate(leftover_unmatched,msms = 'peaks', retention_time='retention_time'):
#     # print("i am in new function")
#     # peaks =[]
#     # leftover_matched = leftover_unmatchedX.drop_duplicates(subset=['key'],keep="first", inplace =False, ignore_index = True)
#     leftover_matched = pd.DataFrame()
#     for key in tqdm(leftover_unmatched['key'].unique()):
#         data_temp = string_search(leftover_unmatched, 'key',key)
#         if len(data_temp)==2:
#             if _deduplicate_duo(data_temp, msms, retention_time)==1:
#                 leftover_matched = pd.concat([leftover_matched, data_temp], axis = 0, ignore_index = True)
#         elif len(data_temp)>2 and (_deduplicate_multi(data_temp,msms, retention_time)) is not None:
#             # peaks.append(_deduplicate_multi(data_temp,msms, retention_time))

#             leftover_matched = pd.concat([leftover_matched, _deduplicate_multi(data_temp,msms, retention_time)], axis = 0, ignore_index = True)
#         elif len(data_temp)==1:
#             print("the length of data_temp is ", len(data_temp))

#             # return(data_temp)


#     # leftover_matched.dropna(subset=['peaks'], inplace = True)
#     return(leftover_matched)


# from toolsets.spectra_operations import entropy_similarity_default
# def _deduplicate_duo(data_temp, msms = 'peaks', retention_time='retention_time'):
#     # print("something")
#     if len(data_temp)!=2:
#         print("this function only handles duplicates with 2")
#         return(0)
#     if entropy_similarity_default(data_temp.iloc[0][msms], data_temp.iloc[1][msms])>=0.75 or abs(data_temp.iloc[0][retention_time]-data_temp.iloc[1][retention_time])<=40/60:
#         # print("i am in if")
#         return(1)
#     else:
#         return(np.NAN)

# def _find_seed_msms(data_temp, msms='peaks', retention_time = 'retention_time'):
#     _count = []
#     for index, row in data_temp.iterrows():
#         entropy_temp = []
#         for index_n, row_n in data_temp.iterrows():
#             entropy_temp.append(entropy_similarity_default(row['peaks'], row_n['peaks']))
#         _count.append(sum(i > 0.75 for i in entropy_temp))
#     if np.max(_count)>1:
#         return(np.argmax(_count))
#         # return(_count)
#     elif data_temp[retention_time].max()-data_temp[retention_time].min() <= 20/60:
#         # print("i am in new if")
#         return(data_temp.reset_index(drop=True)['ms1_precursor_intensity'].idxmax())
#     else:
#         return(-1)




    # return(return_df)

# def _deduplicate_multi(data_temp, msms = 'peaks', retention_time='retention_time'):
#     seed = _find_seed(data_temp, msms = msms)
#     if seed == -1:
#         return(None)
#     picked_spectra = pd.DataFrame()
#     for i in range(0, len(data_temp)):
#         if i != seed:
#             if entropy_similarity_default(data_temp.iloc[seed][msms], data_temp.iloc[i][msms])>=0.75 or abs(data_temp.iloc[seed][retention_time]-data_temp.iloc[i][retention_time])<=6/60:
#                 # print("i am in if")
#                 picked_spectra = picked_spectra.append(data_temp.iloc[i], ignore_index=True)
#                 # print(i)
#             # print("something")
#     picked_spectra = picked_spectra.append(data_temp.iloc[seed], ignore_index=True)
#     return((picked_spectra))

