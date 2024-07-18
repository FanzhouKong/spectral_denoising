import numpy as np
from toolsets.search import quick_search_values, string_search
import pandas as pd
from tqdm import tqdm
import os
from toolsets.file_io import prepare_sample_list, get_file_list
def get_alignemnt(normalized_peak_list_dir, istd_info, mode = 'auto', save_path = None):
    file_list =get_file_list(normalized_peak_list_dir, '.csv', with_tail=True)
    # file_list_no_tail =get_file_list(normalized_peak_list_dir, '.csv', with_tail=False)
    fractions, qcs, blks = prepare_sample_list(file_list)
    print('mode = ', mode)
    if len(qcs)>0:
        print('there is qc, will construct master mz-rt list based on features present in all qcs')
    else:
        print('there is no qc, will construct master mz-rt list from everyone')
    print('starting to construct master mz-rt list')
    if mode == 'auto':
        if len(qcs)>0:
            master_list = initilize_pmz_rt_list(qc_paths=[os.path.join(normalized_peak_list_dir, x) for x in qcs ], mode = 'exclusive')
        else:
            master_list = initilize_pmz_rt_list(qc_paths=[os.path.join(normalized_peak_list_dir, x) for x in fractions ], mode = 'inclusive')
    else:
        if mode =='inclusive':
            if len(qcs)>0:
                master_list = initilize_pmz_rt_list(qc_paths=[os.path.join(normalized_peak_list_dir, x) for x in qcs ], mode = 'inclusive')
            else:
                master_list = initilize_pmz_rt_list(qc_paths=[os.path.join(normalized_peak_list_dir, x) for x in fractions ], mode = 'inclusive')

        elif mode =='exclusive':
            if len(qcs)>0:
                master_list = initilize_pmz_rt_list(qc_paths=[os.path.join(normalized_peak_list_dir, x) for x in qcs ], mode = 'exclusive')
            else:
                print('you dont have any qc, this mode is not supported!')
                return()
    print(f"master mz-rt list done, there is {str(len(master_list))} features in master list!")
    print('aligning features')

    alignment = align(file_paths=[os.path.join(normalized_peak_list_dir, x) for x in file_list ], master_list_input=master_list)
    alignment = string_search(alignment, 'isotopic_state', 'M + 0')
    alignment.drop(['isotopic_state'], axis = 1, inplace = True)
    # alignment_istd = find_istd(alignment, istd_info)

    # if len(blks)>0 and len(qcs)>0:
    #     print('blank files found, cleaning features based on blank')
    #     alignment_istd = filter_with_blank(alignment_result =  alignment_istd, qc_labes=qcs, blank_labels= blks)
    #     print(f"feature cleaning done, there is {str(len(alignment_istd))} features in alignment!")
    # else:
    #     print('no blk or no qc found, no post processing done')
    print('extracting spectra from peak_list')
    alignment_istd = alignment.copy()
    alignment_with_spec = extract_spec(alignment_istd, normalized_peak_dir=normalized_peak_list_dir)
    max_peak_height = []
    for index, row in alignment_with_spec.iterrows():
        max_peak_height.append(row[row['reference_file']])
    alignment_with_spec.insert(2, 'max_peak_height', max_peak_height)
    # alignment_istd.drop(['isotopic_state'], axis = 1, inplace = True)
    # alignment_istd.to_csv(os.path.join(dirs[0], 'alignment_result_istd.csv'), index = False)
    if save_path is not None:
        alignment_with_spec.to_csv(os.path.join(save_path, 'alignment_result.csv'), index = False)
    return(alignment_with_spec)
def align(file_paths, master_list_input):
    master_list = master_list_input.copy()
    master_list_mz_sorted = master_list.sort_values(by = 'pmz', ascending=True)
    master_list_rt_sorted = master_list.sort_values(by = 'rt', ascending=True)
    isotopic_state = np.repeat('unknown', len(master_list))
    for file in tqdm(file_paths):
        file_temp = pd.read_csv(file)
        intensity_temp = np.repeat(0, len(master_list))
        for index, row in file_temp.iterrows():
            rt_matched = quick_search_values(master_list_rt_sorted, 'rt', row['RT_adjusted']-2/60, row['RT_adjusted']+2/60)
            mz_matched = quick_search_values(master_list_mz_sorted, 'pmz', row['Precursor m/z']-0.005, row['Precursor m/z']+0.005)
            lst = list(set(rt_matched.index) & set(mz_matched.index))
            if len(lst)>0:# there is match
                intensity_temp[lst]=intensity_temp[lst]+row['Height']
                if row['Isotope'] =='M + 0':
                    for i in lst:
                        if isotopic_state[i]!='M + 0':
                            isotopic_state[i]='M + 0'


            master_list[os.path.basename(file).split('.')[0]]=intensity_temp
            # feature_temp.drop_duplicates(subset = ['PeakID'])
            # if len(feature_temp)>0:
            #     feature_temp.sort_values(by = ['Height'], inplace = True)
            #     intensity_temp.append(feature_temp.iloc[0]['Height'])
            #     if 'M + 0' in feature_temp['Isotope'].unique() and isotopic_state[index]!= 'M + 0':
            #         # update isotopic state
            #         isotopic_state[index]= 'M + 0'
            # else:
            #     intensity_temp.append(0)

    master_list['isotopic_state']=isotopic_state
    return(master_list)



# def align_legacy(file_paths, master_list_input):
#     master_list = master_list_input.copy()
#     isotopic_state = np.repeat('unknown', len(master_list))
#     for file in tqdm(file_paths):
#         file_temp = pd.read_csv(file)
#         intensity_temp = []
#         # suspicious_df = pd.DataFrame()
#         for index, row in master_list.iterrows():
#             feature_temp = find_feature(file_temp, mz = row['pmz'], rt =row['rt'])
#
#             # feature_temp.drop_duplicates(subset = ['PeakID'])
#             if len(feature_temp)>0:
#                 feature_temp.sort_values(by = ['Height'], inplace = True)
#                 intensity_temp.append(feature_temp.iloc[0]['Height'])
#                 if 'M + 0' in feature_temp['Isotope'].unique() and isotopic_state[index]!= 'M + 0':
#                     # update isotopic state
#                     isotopic_state[index]= 'M + 0'
#             else:
#                 intensity_temp.append(0)
#         master_list[os.path.basename(file).split('.')[0]]=intensity_temp
#     master_list['isotopic_state']=isotopic_state
#     return(master_list)

    # pass
# file_paths = [os.path.join(working_dir, x+'.csv') for x in file_list]
# alignment_result= pd.DataFrame(zip(pmz_list, rt_list), columns=['pmz', 'rt'])
# for file in tqdm(file_paths):
#     file_temp = pd.read_csv(file)
#     intensity_temp = []
#     # pmz_temp = []
#     for i in (range(0, len(pmz_list))):
#         feature_temp = find_feature(file_temp, mz = pmz_list[i], rt =rt_list[i])
#         if len(feature_temp)>0:
#             intensity_temp.append(feature_temp['Height'].sum())
#             # pmz_temp.append()
#             # break
#         else:
#             intensity_temp.append(0)
#     alignment_result[os.path.basename(file).split('.')[0]]=intensity_temp


def initilize_pmz_rt_list(qc_paths,
                          mode='inclusive'
                          ):

    # print('i am in current!')
    if mode == 'inclusive':
        seed_idx = determine_seed_idx(qc_paths)
        qc_seed =  pd.read_csv(qc_paths[seed_idx])

        qc_seed = drop_duplicate_features(qc_seed)
        qc_rest = qc_paths[0:seed_idx]+(qc_paths[seed_idx+1:])
        pmz_list = qc_seed['Precursor m/z'].tolist()
        rt_list = qc_seed['RT_adjusted'].tolist()
        # iso_state = qc_seed['Isotope'].tolist()
        master_list = pd.DataFrame(zip(pmz_list, rt_list), columns=['pmz', 'rt'])
        # return(master_list)
        # print('length_of_master list: ', str(len(master_list)))
        # return(qc_seed)
        if len(qc_rest)>0:
            for qc in tqdm(qc_rest):
                qc_temp = pd.read_csv(qc)
                qc_temp = drop_duplicate_features(qc_temp)
                # break
                for index, row in (qc_temp.iterrows()):
                    feature_temp = find_feature(master_list, mz = row['Precursor m/z'],
                                                rt = row['RT_adjusted'], mz_column='pmz',
                                                rt_column= 'rt'
                                                )
                    if len(feature_temp)<1:
                        pmz_list.append(row['Precursor m/z'])
                        rt_list.append(row['RT_adjusted'])

                master_list = pd.DataFrame(zip(pmz_list, rt_list), columns=['pmz', 'rt'])
                # print('the current length of the master list is: ', str(len(master_list)))
        # print('length_of updated master list: ', str(len(master_list)))

    elif mode == 'exclusive':
        seed_idx = determine_seed_idx(qc_paths, mode = 'exclusive')
        qc_seed =  pd.read_csv(qc_paths[seed_idx])
        qc_seed = drop_duplicate_features(qc_seed)
        # return(qc_seed)
        qc_rest = qc_paths[0:seed_idx]+(qc_paths[seed_idx+1:])
        pmz_list = qc_seed['Precursor m/z'].tolist()
        rt_list = qc_seed['RT_adjusted'].tolist()

        master_list = pd.DataFrame(zip(pmz_list, rt_list), columns=['pmz', 'rt'])
        # return(master_list)
        # print('length_of_master list: ', str(len(master_list)))
        if len(qc_rest)>0:
            for qc in tqdm(qc_rest):
                qc_temp = pd.read_csv(qc)
                qc_temp = drop_duplicate_features(qc_temp)
                drop_idx = []
                for index, row in (master_list.iterrows()):
                    # score.append(distance_score())
                    feature_temp = find_feature(qc_temp, mz = row['pmz'],
                                                rt = row['rt'], mz_column='Precursor m/z',
                                                rt_column= 'RT_adjusted'
                                                )
                    if len(feature_temp)<1:
                        drop_idx.append(index)
                master_list.drop(index = drop_idx, axis = 0, inplace=True)
        # print('length_of updated master list: ', str(len(master_list)))

    else:
        print('wrong mode information is passed, accepted values: inclusive/exclusive')
        return()
    master_list.reset_index(inplace=True, drop = True)
    return master_list
def drop_duplicate_features(peak_listt):
    peak_list = peak_listt.copy()
    peak_list.sort_values(by = ['Height'], ascending=False, inplace=True)
    peak_list['key']=peak_list['Isotope']+peak_list['MS1 isotopes']
    # peak_list['iso_num']= peak_list['Isotope'].str[-1].astype(float)
    # peak_list_return = pd.DataFrame()
    # for key in peak_list['key'].unique():
    #     temp = string_search(peak_list, 'key', key, reset_index=False)
    #     if len(temp)>1:
    #         temp.sort_values(by = 'iso_num', ascending=True, inplace=True)
    peak_list.drop_duplicates(subset = ['key'], inplace = True)
    peak_list.drop(columns=['key'], inplace=True)
    peak_list.sort_values(by = ['RT_adjusted'], inplace=True)
    peak_list.reset_index(inplace=True, drop = True)

    return(peak_list)
def determine_seed_idx(qc_paths, mode = 'inclusive'):
    feature_num = []
    for qc in qc_paths:
        temp= pd.read_csv(qc)
        feature_num.append(len(temp))
    # print(feature_num)
    if mode == 'inclusive':
        return(np.argmax(feature_num))
    else:
        return(np.argmin(feature_num))
def find_istd(alignmentt, istd_info, mz_column = 'pmz', rt_column = 'rt'):
    # print('i am in new')
    alignment = alignmentt.copy()
    istd_idx = []
    for index, row in istd_info.iterrows():
        feature_mz_search = quick_search_values(alignment, mz_column, row['Precursor m/z']-0.005, row['Precursor m/z']+0.005, ifsorted=False)
        # print('i passed first feature search')
        feature_mzrt_search = quick_search_values(feature_mz_search, rt_column, row['RT_suggested']-10/60, row['RT_suggested']+10/60, ifsorted=False)
        # feature_mzrt_search.sort_values(by = '')
        if len(feature_mzrt_search)>1:
            sum_non_zero = []
            for index, row in feature_mzrt_search.iterrows():
                sum_non_zero.append(sum(row[feature_mzrt_search.columns[3:-1]]!=0))
            feature_mzrt_search['sum']=sum_non_zero
            feature_mzrt_search.sort_values(by = 'sum', ascending=False, inplace=True)
            # return(feature_mzrt_search)
        istd_idx.append(feature_mzrt_search.index[0])
        # return(feature_mzrt_search)
    # return istd_idx
    annotation = np.repeat('unknown', len(alignment))
    feature_type = np.repeat('compound', len(alignment))
    feature_type[istd_idx]=np.repeat('istd', len(istd_idx))
    annotation[istd_idx]=istd_info['compound_name']

    alignment.insert(2, 'annotation', annotation)
    alignment.insert(3, 'feature_type', feature_type)
    # alignment['feature_type'].loc[istd_idx]='istd'
    # alignment['annotation'].loc[istd_idx]=

    return (alignment)
import toolsets.spectra_operations as so
def extract_spec(alignment, normalized_peak_dir):
    file_list_no_tail=get_file_list(normalized_peak_dir, '.csv', with_tail=False)
    alignment.insert(2, 'reference_file', alignment[file_list_no_tail].idxmax(axis=1))
    alignment_with_spec = pd.DataFrame()
    for mix in tqdm(alignment['reference_file'].unique()):
        reference_ms1 = []
        reference_msms = []
        alignment_mix = string_search(alignment, 'reference_file', mix)
        feature_list_temp = pd.read_csv(os.path.join(normalized_peak_dir, mix+'.csv'))
        feature_mz_sorted = feature_list_temp.sort_values(by = 'Precursor m/z', ascending=True)
        feature_rt_sorted = feature_list_temp.sort_values(by = 'RT_adjusted', ascending=True)
        for index, row in alignment_mix.iterrows():
            feature_temp = find_feautre_fast(feature_list_temp,feature_mz_sorted,feature_rt_sorted, mz = row['pmz'], rt = row['rt'])
            feature_temp.sort_values(by = 'Height', ascending=False, inplace=True)
            if feature_temp.iloc[0]['MS1 isotopes']==feature_temp.iloc[0]['MS1 isotopes']:
                reference_ms1.append(so.convert_msdial_to_string(feature_temp.iloc[0]['MS1 isotopes'] ))
            else:
                reference_ms1.append(np.NAN)
            if feature_temp.iloc[0]['MSMS spectrum']==feature_temp.iloc[0]['MSMS spectrum']:
                reference_msms.append(so.convert_msdial_to_string(feature_temp.iloc[0]['MSMS spectrum'] ))
            else:
                reference_msms.append(np.NAN)
        alignment_mix.insert(2, 'ms1', reference_ms1)
        alignment_mix.insert(2, 'msms', reference_msms)
        alignment_with_spec = pd.concat([alignment_with_spec, alignment_mix], ignore_index = True)
    return(alignment_with_spec)
        # break
def find_feautre_fast(feature_table, feature_mz_sorted,feature_rt_sorted, mz, rt, mz_column = 'Precursor m/z', rt_column = 'RT_adjusted',  rt_offset = 2):
    step = 0.005
    # feature_mz_sorted = feature_table.sort_values(by = mz_column, ascending=True)
    # feature_rt_sorted = feature_table.sort_values(by = rt_column, ascending=True)
    feature_mz_matched = quick_search_values(feature_mz_sorted, mz_column, mz-step, mz+step, ifsorted=True)
    feature_rt_matched = quick_search_values(feature_rt_sorted, rt_column, rt-rt_offset/60, rt+rt_offset/60, ifsorted=True)
    lst = list(set(feature_mz_matched.index) & set(feature_rt_matched.index))
    return(feature_table.loc[lst])
    # if len(lst)>0:# there is match
    #     return()
def find_feature(feature_table, mz, rt, mz_column = 'Precursor m/z', rt_column = 'RT_adjusted', intensity_column =None, rt_offset = 2):
    # print('i am nin new')
    mz_step = 0.003
    # print(mz_step)
    feature_mz_search = quick_search_values(feature_table, mz_column, mz-mz_step, mz+mz_step, ifsorted=False)
    # return(feature_mz_search)
    feature_mzrt_search = quick_search_values(feature_mz_search, rt_column, rt-rt_offset/60, rt+rt_offset/60, ifsorted=False)
    if intensity_column != None:
        feature_mzrt_search.sort_values(by = intensity_column, inplace=True, ascending=False)
    # print(feature_mzrt_search)
    return (feature_mzrt_search)
def clean_bad_features(alignment_result):
    alignment_result_refined = alignment_result.copy()
    bad_idx = []
    for index, row in alignment_result.iterrows():
        if np.sum(row[4:]) == np.max(row[4:]):
            bad_idx.append(index)
    alignment_result_refined.drop(bad_idx, inplace = True)
    alignment_result_refined.reset_index(inplace = True, drop = True)
    return(alignment_result_refined)
from toolsets.file_io import get_list_idx
def filter_with_blank(alignment_result, qc_labes, blank_labels, threshold = 1, mode = 'median'):
    # to make this consistent, all qc_labels/blk_labels coming in there still has tail .csv
    # so the first step is to remove them
    qc_labes = [[x[0:-4] for x in qc_labes]]
    blank_labels = [[x[0:-4] for x in blank_labels]]
    alignment = alignment_result.copy()

    alignment_result_qcs = alignment[qc_labes]
    alignment_result_blk = alignment[blank_labels]
    bad_feature_idx = []
    for index, row in alignment.iterrows():
        if mode == 'median':
            toggle = alignment_result_qcs.loc[index].median()*threshold/100<alignment_result_blk.loc[index].max()
        elif mode =='max':
            toggle = alignment_result_qcs.loc[index].max()*threshold/100<alignment_result_blk.loc[index].max()
        if row['feature_type']!='istd' and toggle == True:
            bad_feature_idx.append(index)
    alignment.drop(bad_feature_idx, inplace = True)
    alignment.drop(blank_labels, axis = 1, inplace = True)
    alignment.reset_index(inplace=True, drop = True)
    return(alignment)
    # check if blank labels is given
    # if blank_lables is None:

    # blank_df = alignment_result[blank_lables]