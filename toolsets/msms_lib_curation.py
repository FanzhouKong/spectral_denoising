
from toolsets.ff_droup import process_mzml, feature_finding
from toolsets.search import quick_search_values, quick_search_sorted, string_search
import os
import pandas as pd
import numpy as np
import toolsets.ff_droup as ff
from tqdm import tqdm
import toolsets.spectra_operations as so
def extract_features(mzml_dir, features_dir, rt_max = 20):
    if os.path.exists(features_dir)==False:
        os.mkdir(features_dir)
    for root, dirs, files in os.walk(mzml_dir):
        for file in tqdm(files):
            if file.endswith('.mzML'):
                try:
                    ms1, ms2 =ff.process_mzml(os.path.join(mzml_dir, file), if_mix=False, rt_max=rt_max)
                    features_temp = ff.feature_finding(ms1, ms2)
                    features_temp.to_csv(os.path.join(features_dir, file.split('.')[0]+'.csv'), index = False)
                except:
                    print('this file is worng: ', file)
def export_msp_simple_msp(data,output_location):
    entry = ''
    for index, row in tqdm(data.iterrows(), total = len(data)):
        entry = entry + 'Name: ' + row['reference_name'] + '\n'
        entry = entry + 'RETENTIONTIME: ' + str(row['rt']) + '\n'
        entry = entry +'Spectrum_type: '+'MS2'+ '\n'
        entry = entry + 'PrecursorMZ: ' + str(row['precursor_mz']) + '\n'
        entry = entry +'Precursor_type: ' + str(row['reference_adduct']) + '\n'
        entry = (entry + 'Comment: ' + 'ms1intensity'+'_'+str(row['peak_apex_intensity'])+ '\n')
        entry = entry + 'Num peaks: ' + str(so.num_peaks(row['peaks'])) + '\n'
        entry = entry + row['peaks']
        entry = entry + '\n'
        entry = entry + '\n'

    #open text file
    text_file = open(output_location, "w",encoding='utf-8')

    #write string to file
    text_file.write(entry)

    #close file
    text_file.close()
from toolsets.feature_alignment import find_feature
def library_curation(std_list, features_dir, adducts = ['[M+H]+','[M+Na]+','[M-H2O+H]+','[M+K]+']):
    matched = pd.DataFrame()
    bad_mix = []
    matched= pd.DataFrame()
    for mix in tqdm(std_list['mix'].unique()):
        try:
            feature_all = pd.read_csv(os.path.join(features_dir, mix+'.csv'))
        except:
            print('this mix is wrong')
            continue
        std_mix = string_search(std_list, 'mix', mix)

        for index, row in std_mix.iterrows():
            if 'rt' in std_list.columns.str.lower() or 'retention time' in std_list.columns.str.lower():
                pass
            else:
                for adduct in adducts:
                    feature_temp = quick_search_values(feature_all, 'precursor_mz', float(row[adduct])-0.003,
                                                           float(row[adduct])+0.003, ifsorted=False)
                    if len(feature_temp)>0:# there is a match
                        feature_temp['rt_offset']=abs(feature_temp['rt_offset'])
                        feature_temp.sort_values(by = 'rt_offset', ascending=True, inplace=True)
                        feature_temp.insert(0, 'reference_name', np.repeat(row['Name'], len(feature_temp)))
                        feature_temp.insert(1, 'reference_precursor_mz', np.repeat(row[adduct], len(feature_temp)))
                        feature_temp.insert(1, 'reference_smiles', np.repeat(row['SMILES'], len(feature_temp)))
                        feature_temp.insert(2,'reference_mix',np.repeat(str(row['BSD_ID']), len(feature_temp)))
                        feature_temp.insert(2, 'reference_adduct', np.repeat(adduct, len(feature_temp)))
                        matched = pd.concat([matched, pd.DataFrame([feature_temp.iloc[0]])], ignore_index = True)
    matched_filtered = pd.DataFrame()
    missing_compound = list(set(std_list['Name']) - set(matched['reference_name']))
    for name in matched['reference_name'].unique():
        compound_temp = string_search(matched, 'reference_name', name)
        compound_temp.sort_values(by = 'peak_apex_intensity', ascending=False, inplace=True)
        seed_rt = compound_temp.iloc[0]['rt']
        compound_filtered = quick_search_values(compound_temp, 'rt', seed_rt-5/60, seed_rt+5/60, ifsorted=False)
        matched_filtered = pd.concat([matched_filtered, compound_filtered], ignore_index=True)
    # bad_mix = [x for x in std_list['mix'] if x not in matched['reference_mix']]
    return(matched_filtered, missing_compound)


from toolsets.API_gets import name_to_smiles
from toolsets.std_list_prep import calculate_precursormz
from toolsets.constants import single_charged_adduct_mass
def complete_std_list(std_listt, adducts = ['[M+H]+','[M+Na]+','[M-H2O+H]+','[M+K]+'], mode = 'pos'):
    std_list = std_listt.copy()
    smiles = []
    mix = []
    # print('i am fetching smiles')
    for index, row in (std_list.iterrows()):
        if row['SMILES'] != row['SMILES']:
            smiles.append(name_to_smiles(row['Name']))
        else:
            smiles.append(row['SMILES'])
    # print('i start calculating adducts')
    std_list['SMILES']=smiles
    std_list.dropna(subset=['SMILES'], inplace=True, ignore_index=True)
    for adduct in adducts:
        pmz_adduct = []
        for index, row in std_list.iterrows():
            pmz_adduct.append(calculate_precursormz(row['SMILES'], adduct))
        std_list[adduct]=pmz_adduct
    if mode == 'pos':
        std_list['mix']=std_list['BSD_ID']+'_P'
    elif mode == 'neg':
        std_list['mix']=std_list['BSD_ID']+'_N'
    else:
        std_list['mix']=std_list['BSD_ID']
    return(std_list)