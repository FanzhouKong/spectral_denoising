import pandas as pd
from toolsets.search import  string_search, quick_search_values, quick_search_sorted
from toolsets.spectra_operations import entropy_identity
import numpy as np
def dereplicate(compound_matched):
    if len(compound_matched)==0:
        return(pd.DataFrame())
    df_return = pd.DataFrame()
    guessed_rt = compound_matched.iloc[np.argmax(compound_matched['ms1_intensity'])]['rt_apex']
    compound_matched.dropna(subset = 'msms', inplace = True)
    if len(compound_matched)==0:
        return(compound_matched)
    comment = []
    for ma in compound_matched['reference_adduct'].unique():
        # comment = []
        current_adduct = string_search(compound_matched, 'reference_adduct', ma)
        rt_matched = quick_search_values(current_adduct, 'rt_apex', guessed_rt-5/60, guessed_rt+5/60)
        if len(rt_matched)>0:
            current_adduct_left_over = current_adduct.drop(rt_matched.index.tolist())
            major = rt_matched[rt_matched['ms1_intensity']==rt_matched['ms1_intensity'].max()]
            minor = rt_matched[rt_matched['ms1_intensity']<rt_matched['ms1_intensity'].max()]
            if len(major)>1:
                major.sort_values(by = 'rt_offset', ascending=True, inplace = True)
                major = major[0:1]
            df_return = pd.concat([df_return, major], ignore_index=True)
            comment.append('Major')
            df_return = pd.concat([df_return, minor], ignore_index=True)
            comment.extend(['Minor']*len(minor))
            seed_msms = major.iloc[0]['msms']
            for i, j in current_adduct_left_over.iterrows():
                entropy_temp =entropy_identity(seed_msms, j['msms'], pmz = major.iloc[0]['precursor_mz'])
                if entropy_temp>0.75:
                    df_return= pd.concat([df_return, pd.DataFrame([j])], ignore_index=True)
                    comment.append('isomer')
    if len(df_return)==0:
        return (df_return)
    df_return.insert(4, 'comment', comment)
    return df_return
def feature_matching(feature_mix, std_list_mix, adducts, unique_identifier = 'CAS'):
    feature_mix.sort_values(by = 'precursor_mz', inplace=True, ascending=True)
    mix_matched = pd.DataFrame()
    for index, row in std_list_mix.iterrows():
        compound_matched = pd.DataFrame()

        for a in adducts:
            adduct_matched = quick_search_sorted(feature_mix, 'precursor_mz', row[a]-0.005, row[a]+0.005)
            if len(adduct_matched)>0:
                # adduct_matched.insert(0, 'reference_mz', row[a])
                adduct_matched.insert(1, 'reference_name', row['name'])
                adduct_matched.insert(1, 'reference_mz', row[a])
                adduct_matched.insert(1,unique_identifier, row[unique_identifier])
                adduct_matched.insert(2, 'reference_adduct', a)

                # adduct_matched.insert(3, 'reference_rt', row['reference_rt'])
                adduct_matched.insert(4, 'reference_smiles', row['smiles'])
                adduct_matched.insert(5, 'reference_formula', row['formula'])
                adduct_matched.insert(6, 'reference_mix', row['mix'])
                adduct_matched.insert(7, 'reference_rt', row['rt_reference'])
                compound_matched  = pd.concat([compound_matched, adduct_matched], ignore_index=True)
        if len(compound_matched)>0:
            compound_matched = dereplicate(compound_matched)
            mix_matched = pd.concat([mix_matched, compound_matched],ignore_index=True)
            # compound_matched = pd.concat([compound_matched, adduct_matched], ignore_index=True)




    if len(mix_matched)>0:
        mix_matched.drop(columns=['eic_center', 'eic_offset'], inplace=True)
    return(mix_matched)