from toolsets.helpers import find_floats
import numpy as np
import toolsets.chem_utils as ch
def complete_std_list(std_list, smiles_col = 'smiles', adducts = ['[M]+', '[M+H]+', '[M+Na]+','[M+NH4]+']):
    std_list = complete_formula(std_list, smiles_col = smiles_col)
    std_list = complete_mono_mass(std_list, smiles_col = smiles_col)
    std_list = complete_adducts(std_list, smiles_col = smiles_col, adducts = adducts)
    return(std_list)
def clean_rt(std_list, colname = 'rt_suggested'):
    df = std_list.copy()
    rt_cleaned = []
    for r in df[colname]:
        if r == r:
            rt_cleaned.append(find_floats(r))
        else:
            rt_cleaned.append(np.NAN)
    df[colname]=rt_cleaned
    return df
def complete_formula(std_list_raw, smiles_col = 'smiles'):
    std_list = std_list_raw.copy()
    formulas = []
    for index, row in std_list.iterrows():
        if ch.is_smiles(row[smiles_col])==True:
            formulas.append(ch.everything_to_formula(row[smiles_col]))
        else:
            formulas.append(np.NAN)
    column_idx = std_list.columns.get_loc(smiles_col)
    std_list.insert(column_idx+1, 'formula',formulas)
    return std_list
def complete_mono_mass(std_list_raw, smiles_col = 'smiles'):
    std_list = std_list_raw.copy()
    mono_mass = []
    for index, row in std_list.iterrows():
        if ch.is_smiles(row[smiles_col])==True:
            mono_mass.append(ch.everything_to_mw(row[smiles_col]))
        else:
            mono_mass.append(np.NAN)
    std_list['mono_mass']=mono_mass
    return std_list
def complete_adducts(std_list_raw, smiles_col = 'smiles', adducts = '[M+H]+'):
    std_list = std_list_raw.copy()
    adducts_dict = {}
    for adduct in adducts:
        adducts_dict[adduct]=[]
    for index, row in std_list.iterrows():
        if ch.is_smiles(row[smiles_col])==True:
            for adduct in adducts:
                adducts_dict[adduct].append(ch.calculate_precursormz(row[smiles_col], adduct=adduct, if_smiles=True))
        else:
            for adduct in adducts:
                adducts_dict[adduct].append(np.NAN)
    for adduct in adducts:
        std_list[adduct]=adducts_dict[adduct]
    return std_list