import pandas as pd
from tqdm import tqdm
import itertools
import numpy as np
import os

from toolsets.search import string_search, quick_search_sorted
from toolsets.spectra_operations import normalize_spectrum, spectral_entropy, clean_spectrum

def precursor_matching_multiple(feature_paths, std_list_paths, library_names):
    if len(feature_paths)!= len(std_list_paths) or len(std_list_paths)!= len(library_names):
        raise RuntimeError
    matched_all = pd.DataFrame()
    for i in range(0, len(feature_paths)):
        feature = pd.read_csv(feature_paths[i])
        std_list = pd.read_csv(std_list_paths[i])
        matched_temp= precursor_matching_dataframe(feature, std_list)
        matched_temp['subdivision']=library_names[i]
        matched_all = pd.concat([matched_all, matched_temp], axis=0, ignore_index=True)
    matched_all['library_id'] = np.arange(len(matched_all))
    return(matched_all)
def precursor_matching_dataframe(features_all, std_list_all, if_clean = True):
    matched = pd.DataFrame()
    for mix in tqdm(std_list_all['mix'].unique()):
        feature_mix = string_search(features_all, 'mix', mix)
        std_list_mix = string_search(std_list_all, 'mix', mix)
        matches_temp = _precursor_matching_mix(std_list_mix, feature_mix)
        matched = pd.concat([matched, matches_temp], axis = 0)
    # return (matched)

    print("i have that much spectra before cleaninig:", str(len(matched)))
    print("i have that much compounds before cleaning:", str(len(matched['reference_inchikey'].unique())))
    if if_clean ==True:
        peaks_cleaned = []
        spectral_entropy_lst = []
        for index, row in tqdm(matched.iterrows(), total=len(matched)):
            cleaned_peak = clean_spectrum(row['peaks'], row['reference_precursor_mz'])
            peaks_cleaned.append(cleaned_peak)
            spectral_entropy_lst.append(spectral_entropy(cleaned_peak))

        matched['peaks']=peaks_cleaned
        matched['spectral_entropy']=spectral_entropy_lst
        matched = matched[matched['spectral_entropy']>0.5]
        matched.reset_index(inplace=True, drop=True)
        matched['library_id'] = np.arange(len(matched))

    return(matched)

def find_adduct_columns(columns, adducts = None):
    from toolsets.constants import single_charged_adduct_mass
    # from toolsets.std_list_prep import all_adduct_pos, all_adduct_neg
    if adducts == None:
        adducts = []
    else:
        adducts = adducts
    for col in columns:
        if col in single_charged_adduct_mass.keys():
            adducts.append(col)
    return(adducts)
def _precursor_matching_mix(std_list_mix, features, step = 0.005):
    features.sort_values(by='precursor_mz', inplace=True, ignore_index=True)
    adducts = find_adduct_columns(std_list_mix)
    matched_mix = pd.DataFrame()
    raw_mix = pd.DataFrame()
    reference_columns = [col for col in std_list_mix.columns if col not in adducts]
    reference_adduct = []
    reference_precursor_mz=[]
    for index, row in std_list_mix.iterrows():
        for adduct in adducts:
            # print(row[adduct])

            raw_matched = quick_search_sorted(features, 'precursor_mz', float(row[adduct])-step, float(row[adduct])+step)

            if len(raw_matched)>0:   #there is at least 1 valid match
                # raw_matched.sort_values(by = 'ms1_intensity', inplace=True, ascending= False)

                matched_mix= matched_mix.append([pd.DataFrame([row[reference_columns]])]*len(raw_matched), ignore_index=True)
                raw_mix= pd.concat([raw_mix, raw_matched],ignore_index=True)
                reference_adduct.extend([adduct]*len(raw_matched))
                reference_precursor_mz.extend([row[adduct]]*len(raw_matched))
        #     break
        # break
    # print(len(matched_mix))
    # return (matched_mix)
    matched_mix.columns = ['reference_'+col_name for col_name in matched_mix.columns]
    matched_mix['reference_adduct']= reference_adduct
    matched_mix['reference_precursor_mz']=reference_precursor_mz
    matched_mix_spec = pd.concat([matched_mix, raw_mix], axis=1)

    if len(matched_mix_spec)>0:
        matched_mix_spec['mz_offset'] = abs(matched_mix_spec['precursor_mz']-matched_mix_spec['reference_precursor_mz'])
        recovery_percent = len(set(matched_mix_spec['reference_inchikey']))/len(set(std_list_mix['inchikey']))*100
        matched_mix_spec['recovery_percent(%)']=recovery_percent
    return(matched_mix_spec)
def _precursor_matching_mix_msdial(std_list_mix, height, step = 0.005):
    height.sort_values(by='Average Mz', inplace=True, ignore_index=True)
    # height.drop(subset = [])
    adducts = find_adduct_columns(std_list_mix)
    matched_mix = pd.DataFrame()
    raw_mix = pd.DataFrame()
    reference_columns = [col for col in std_list_mix.columns if col not in adducts]
    reference_adduct = []
    reference_precursor_mz=[]
    for index, row in std_list_mix.iterrows():
        for adduct in adducts:
            raw_matched = quick_search_sorted(height, 'Average Mz', float(row[adduct]), step = step)
            if len(raw_matched)>0:   #there is at least 1 valid match
                # raw_matched.sort_values(by = 'ms1_intensity', inplace=True, ascending= False)
                matched_mix= pd.concat([matched_mix, pd.DataFrame([row[reference_columns]])], ignore_index=True)
                raw_mix= pd.concat([raw_mix, pd.DataFrame([raw_matched.iloc[0]])],ignore_index=True)
                reference_adduct.append(adduct)
                reference_precursor_mz.append(row[adduct])
        #     break
        # break
    matched_mix.columns = ['reference_'+col_name for col_name in matched_mix.columns]
    matched_mix['reference_adduct']= reference_adduct
    matched_mix['reference_precursor_mz']=reference_precursor_mz
    matched_mix_spec = pd.concat([matched_mix, raw_mix], axis=1)
    # matched_mix_spec['mz_offset'] = abs(matched_mix_spec['precursor_mz']-matched_mix_spec['reference_precursor_mz'])
    # if len(matched_mix_spec)>0:
    #     recovery_percent = len(set(matched_mix_spec['reference_inchikey']))/len(set(std_list_mix['inchikey']))*100
    #     matched_mix_spec['recovery_percent(%)']=recovery_percent
    return(matched_mix_spec)

def get_bootstrap_std():
    pass
# def EIC(pmz, ms1):

# def precursor_matching_global(std_list, feature_dir, mass_error = 0.005, ppm = False, comment = None, ifrt = False):
#     adducts = find_adduct_columns(std_list.columns)
#     matched = pd.DataFrame()
#     for mix in (std_list['mix'].unique()):
#         feature_filename = find_files(feature_dir, mix+'.csv')
#         if len(feature_filename)==1:
#             feature = read_in_alphapept(os.path.join(feature_dir, feature_filename[0]))
#             matched_temp = precursor_matching_alphapept(feature, std_list, mix_label=mix, adducts = adducts,step=mass_error, ppm=ppm, ifrt= ifrt)
#             matched = pd.concat([matched, matched_temp], axis=0)
#         else:
#             print("there is either multiple or none feature file for this %s" % str(mix))


#     matched['key']= matched['reference_inchikey']+matched['reference_adduct']
#     # return(matched)
#     spectrum_entropy = []
#     for index, row in matched.iterrows():
#         spectrum_entropy.append(spectral_entropy(row['peaks']))
#     matched['spectrum_entropy']=spectrum_entropy
#     matched = num_search(matched, 'spectrum_entropy', 0.5, direction = '>', inclusion= True)

#     # return(matched)
#     print("cleaning spectrum...")
#     cleaned_peaks = []
#     for index, row in tqdm(matched.iterrows(), total = len(matched)):
#         cleaned_peaks.append(clean_spectrum(msms = row['peaks'],
#                                             max_mz = row['precursor_mz'],
#                                             tolerance = 0.02,
#                                             ifppm = False,
#                                             noise_level= 0.000))
#     matched['peaks']=cleaned_peaks
#     matched.dropna(subset = ['peaks'], inplace = True)
#     matched.reset_index(inplace=True, drop=True)
#     if comment != None:
#         matched['comment']=comment
#     else:
#         matched['comment']=feature_dir

#     # matched['raw_id'] = np.arange(len(matched))
#     spectrum_entropy = []
#     for index, row in matched.iterrows():
#         spectrum_entropy.append(spectral_entropy(row['peaks']))
#     matched = num_search(matched, 'spectrum_entropy', 0.5, direction = '>', inclusion= True)
#     matched.dropna(subset = ['peaks'], inplace = True)
#     matched.reset_index(inplace=True, drop=True)
#     return matched

# def precursor_matching_alphapept(features, std_list, mix_label, adducts = ['[M+H]+'], step = 0.01, ppm = False,
#     information_column = ['name','inchikey','mix','smiles','formula','mono_mass','formal_charges'], ifrt = False

#     ):
#     if ifrt == True:
#         information_column.append('rt')
#     # features = read_in_alphapept(features_path)
#     # std_list = string_search(std_list, "Mix label", mix_label)
#     std_list = string_search(std_list, "mix", mix_label)
#     picked = pd.DataFrame()
#     reference_name = []

#     for index, row in std_list.iterrows():
#         # print("the input parameters are step is %s, ppm is %s" %(str(step), str(ppm)))
#         if ppm == True:
#             if row['mono_mass']<=400:
#                 step_use = 0.004
#             else:
#                 step_use = row['mono_mass']*step/1E6
#         else:
#             step_use = step

#         # print("the monoiso is ", row['Monoisotopic mass'])
#         # print("the step is", (step_use))
#         for adduct in adducts:
#             temp = num_search(features, "precursor_mz", row[adduct], 'between', step =step_use, inclusion=True)
#             if len(temp)>0:
#                 for col in information_column:
#                     temp_col = "reference_"+col
#                     temp[temp_col]=row[col]
#                 temp['reference_adduct']=adduct
#                 temp['reference_precursor_mz']=row[adduct]
#                 picked = pd.concat([picked, temp], axis=0)
#     # picked = num_search(picked, "ms1_intensity_ratio", lowest_ms1_ratio, direction = ">", inclusion = True)
#     return(picked)


# def GT_precursor_matching(std_list, enzyme_list,feature_dir,  adducts = ['[M+H]+'], step = 0.01, ppm = False):
#     matched = pd.DataFrame()
#     print("i am in new method")
#     for mix in tqdm(std_list['Mix'].unique()):
#         for enzyme in (enzyme_list['List of enzymes'].unique()):
#             filename =enzyme+"_"+"Mix"+"_"+str(mix)
#             feature_path = os.path.join(feature_dir, filename+".csv")
#             if os.path.exists(feature_path) ==True:
#                 try:
#                     features = read_in_alphapept(feature_path)
#                     matched_temp = precursor_matching_alphapept(features, std_list, mix_label=mix, adducts = adducts,step=step, ppm=ppm)
#                     # return(matched_temp)
#                     matched_temp['reference_mix_enzyme'] = filename
#                     matched = pd.concat([matched, matched_temp], axis=0)
#                 except:
#                     print("the problematic mix is mix %s and enzyme %s" %(str(mix), str(enzyme)))
#                     # return(features)

#             else:
#                 # print("the feature set for mix %s does not exist"% (mix))
#                 pass
#     matched.reset_index(inplace=True, drop=True)
#     # matched['key']= matched['reference_InChIKey']+matched['reference_adduct']
#     return(matched)
# def precursor_matching(df, sample_list, mzmmu, adducts, comments, CE = None,threshold = 0.3, ifppm = False, iffill = True):
#     # print("i am updated!")
#     # mz_mmu = mzmmu
#     data_raw = pd.DataFrame(columns = ['NAME','key','PRECURSORMZ','InChIKey','Formula',
#                        'ExactMass','Adduct','Spectrum_type','RETENTIONTIME','Average_mz',
#                        'Comment', 'Alignment_ID','ms1','msms','Collision_energy','intensity','mix_label'])
#     if iffill == True:
#         ms2df = df.loc[(df['MS/MS assigned']) & (df['Fill %'] < threshold)]
#     else:
#         ms2df = df.loc[(df['MS/MS assigned'])]
#     # parameters
#     # this is already matching the spectrum in the same data file
#     spectrum_type = 'MS2'
#     # adducts = ['[M+H]+', '[M+NH4]+', '[M+Na]+']
#     # comments = ' EAD spectrum'
#     for mix in tqdm(ms2df.columns[32:-2]):
#         for index, df_row in ms2df[ms2df['Spectrum reference file name'] == mix].iterrows():
#             msms = df_row['MS/MS spectrum'].replace(' ', '\n').replace(':', '\t')
#             ms1 = df_row['MS1 isotopic spectrum'].replace(' ', '\n').replace(':', '\t')
#             nlines = msms.count('\n')+1
#             if ifppm:
#                 if df_row['Average Mz']<400:
#                     mz_lower= df_row['Average Mz'] - 0.004
#                     mz_upper = df_row['Average Mz'] + 0.004
#                 else:
#                     mz_lower = df_row['Average Mz'] - df_row['Average Mz']*mzmmu/(1E6)
#                     mz_upper = df_row['Average Mz'] + df_row['Average Mz']*mzmmu/(1E6)
#             else:
#                 mz_lower= df_row['Average Mz'] - mzmmu
#                 mz_upper = df_row['Average Mz'] + mzmmu

#             for index, ls_row in sample_list[sample_list['Mix label'] == mix].iterrows():
#                 for adduct in adducts:
#                     if ls_row[adduct] > mz_lower  and ls_row[adduct] < mz_upper:
#                         if adduct in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
#                             ion_mode = 'P'
#                         else:
#                             ion_mode = 'N'
#                     if ls_row[adduct] > mz_lower  and ls_row[adduct] < mz_upper:
#                         data_raw = data_raw.append({
#                                                         'NAME':ls_row['Name'],
#                                                        'key':str(ls_row['InChIKey'])+adduct,
#                                                         'Ion_mode':ion_mode,
#                                                        'PRECURSORMZ':(float(ls_row[adduct])),
#                                                        'InChIKey':ls_row['InChIKey'],
#                                                         'Formula':ls_row['Formula'],
#                                                         'ExactMass':(float(ls_row['Exact Mass'])),
#                                                         'Adduct':adduct,
#                                                         'Collision_energy':CE,
#                                                         'Spectrum_type':spectrum_type,
#                                                         'RETENTIONTIME':(float(df_row['Average Rt(min)'])),
#                                                         'Average_mz':(float(df_row['Average Mz'])),
#                                                         'Comment':str(df_row['Alignment ID']) + comments + ' intensity ' + str(df_row[mix]),
#                                                         'Alignment_ID':(df_row['Alignment ID']),
#                                                         # 'Num_Peaks':int(nlines),
#                                                         'ms1':ms1,
#                                                         'msms':msms,
#                                                         'intensity':(df_row[mix]),
#                                                         'mix_label':mix,
#                                                        # 'count':str('na')
#                             }, ignore_index=True)
#     # occ = data_raw.key.value_counts(normalize=False, sort = False)


#     # data_raw_count = pd.DataFrame()
#     # for i in range(len(occ)):
#     #     temp_df = data_raw.loc[data_raw.key == occ.index[i],:]
#     #     temp_df['count'] = temp_df['count'].replace(['na'], int(occ[i]))
#     #     data_raw_count=pd.concat([data_raw_count, temp_df], ignore_index = True, axis = 0)
#     # data_raw_count['row_num'] = np.arange(len(data_raw_count))
#     data_raw.round({'PRECURSORMZ':6})
#     # data_raw['PRECURSORMZ']=pd.to_numeric(data_raw['PRECURSORMZ'])
#     # data_raw['ExactMass']=pd.to_numeric(data_raw['ExactMass'])
#     # data_raw['Average_mz']=pd.to_numeric(data_raw['Average_mz'])
#     # data_raw['intensity']=pd.to_numeric(data_raw['intensity'])
#     # data_raw['Num_Peaks']=pd.to_numeric(data_raw['Num_Peaks'])
#     # data_raw['RETENTIONTIME']=pd.to_numeric(data_raw['RETENTIONTIME'])
#     data_raw.reset_index(inplace=True, drop=True)
#     return(data_raw)

# def precursor_matching_rt(df, sample_list, mzmmu, adducts, comments, CE = None,threshold = 0.3, ppm = False):
#     # print("i am updated!")
#     # mz_mmu = mzmmu
#     data_raw = pd.DataFrame(columns = ['NAME','key','PRECURSORMZ','InChIKey','Formula',
#                        'ExactMass','Adduct','Spectrum_type','RETENTIONTIME','Average_mz',
#                        'Comment', 'Alignment_ID','ms1','msms','Collision_energy','intensity','mix_label'])
#     ms2df = df.loc[(df['MS/MS assigned']) & (df['Fill %'] < threshold)]
#     # parameters
#     # this is already matching the spectrum in the same data file
#     spectrum_type = 'MS2'
#     # adducts = ['[M+H]+', '[M+NH4]+', '[M+Na]+']
#     # comments = ' EAD spectrum'
#     for mix in tqdm(ms2df.columns[32:-2]):
#         for index, df_row in ms2df[ms2df['Spectrum reference file name'] == mix].iterrows():
#             msms = df_row['MS/MS spectrum'].replace(' ', '\n').replace(':', '\t')
#             ms1 = df_row['MS1 isotopic spectrum'].replace(' ', '\n').replace(':', '\t')
#             nlines = msms.count('\n')+1
#             if ppm:
#                 if df_row['Average Mz']<400:
#                     mz_lower= df_row['Average Mz'] - 0.004
#                     mz_upper = df_row['Average Mz'] + 0.004
#                 else:
#                     mz_lower = df_row['Average Mz'] - df_row['Average Mz']*mzmmu/(1E6)
#                     mz_upper = df_row['Average Mz'] + df_row['Average Mz']*mzmmu/(1E6)
#             else:
#                 mz_lower= df_row['Average Mz'] - mzmmu
#                 mz_upper = df_row['Average Mz'] + mzmmu
#             rt = df_row['Average Rt(min)']
#             for index, ls_row in sample_list[sample_list['Mix label'] == mix].iterrows():
#                 for adduct in adducts:
#                     if ls_row[adduct] > mz_lower  and ls_row[adduct] < mz_upper:
#                         if adduct in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
#                             ion_mode = 'P'
#                         else:
#                             ion_mode = 'N'
#                     if ls_row[adduct] > mz_lower  and ls_row[adduct] < mz_upper:
#                         if abs(ls_row['RT [min]']-rt)<0.5:
#                             data_raw = data_raw.append({
#                                                         'NAME':ls_row['Name'],
#                                                        'key':ls_row['InChIKey']+adduct,
#                                                         'Ion_mode':ion_mode,
#                                                        'PRECURSORMZ':(float(ls_row[adduct])),
#                                                        'InChIKey':ls_row['InChIKey'],
#                                                         'Formula':ls_row['Formula'],
#                                                         'ExactMass':(float(ls_row['Exact Mass'])),
#                                                         'Adduct':adduct,
#                                                         'Collision_energy':CE,
#                                                         'Spectrum_type':spectrum_type,
#                                                         'RETENTIONTIME':(float(df_row['Average Rt(min)'])),
#                                                         'Average_mz':(float(df_row['Average Mz'])),
#                                                         'Comment':str(df_row['Alignment ID']) + comments + ' intensity ' + str(df_row[mix]),
#                                                         'Alignment_ID':(df_row['Alignment ID']),
#                                                         # 'Num_Peaks':int(nlines),
#                                                         'ms1':ms1,
#                                                         'msms':msms,
#                                                         'intensity':(df_row[mix]),
#                                                         'mix_label':mix,
#                                                        # 'count':str('na')
#                             }, ignore_index=True)

#     # occ = data_raw.key.value_counts(normalize=False, sort = False)


#     # data_raw_count = pd.DataFrame()
#     # for i in range(len(occ)):
#     #     temp_df = data_raw.loc[data_raw.key == occ.index[i],:]
#     #     temp_df['count'] = temp_df['count'].replace(['na'], int(occ[i]))
#     #     data_raw_count=pd.concat([data_raw_count, temp_df], ignore_index = True, axis = 0)
#     # data_raw_count['row_num'] = np.arange(len(data_raw_count))
#     data_raw.round({'PRECURSORMZ':6})
#     # data_raw['PRECURSORMZ']=pd.to_numeric(data_raw['PRECURSORMZ'])
#     # data_raw['ExactMass']=pd.to_numeric(data_raw['ExactMass'])
#     # data_raw['Average_mz']=pd.to_numeric(data_raw['Average_mz'])
#     # data_raw['intensity']=pd.to_numeric(data_raw['intensity'])
#     # data_raw['Num_Peaks']=pd.to_numeric(data_raw['Num_Peaks'])
#     # data_raw['RETENTIONTIME']=pd.to_numeric(data_raw['RETENTIONTIME'])
#     return(data_raw)
# print("I am precursor matching!")