import os
import sys


import numpy as np
import re
import chemparse
from rdkit.Chem.MolStandardize import rdMolStandardize
from molmass import Formula
import itertools
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog('rdApp.*')
from toolsets.chem_utils import parse_ind, calculate_precursormz, break_adduct, break_adduct
import toolsets.spectra_operations as so
from toolsets.search import quick_search_sorted, string_search, quick_search_values
import pandas as pd
def denoising_search(msms, pmz, library, mass_error = 0.01, pmz_col = 'precursor_mz', smiles_col = 'smiles', adduct_col = 'adduct', spectrum_col = 'spectrum'):
    raw_msms = (msms)
    pmz_candidate = quick_search_values(library, pmz_col, pmz-mass_error, pmz+mass_error)
    charge = [break_adduct(x)[2] for x in pmz_candidate[adduct_col]]
    processed_formulas = [prep_formula(x[smiles_col], x[adduct_col]) for i, x in pmz_candidate.iterrows()]
    pmz_candidate['processed_formulas']=processed_formulas
    pmz_candidate['charge']=charge
    entropy_denoised = []
    entropy_raw = []
    msms_denoised_all = []
    pmz_candidate_curated = pd.DataFrame()
    for f in pmz_candidate['processed_formulas'].unique():

        if f ==f:
            pmz_f_temp = string_search(pmz_candidate, 'processed_formulas', f)
            element_dict, element_count, element_mass, all_possible_mass, all_possible_candidate_formula = get_all_subsets(f)

            msms_denoised = denoising_with_subset(raw_msms, pmz, element_dict, all_possible_mass, all_possible_candidate_formula)
            msms_denoised_temp = [msms_denoised]*len(pmz_f_temp)
            entropy_raw_temp = [so.entropy_identity(msms, m[spectrum_col], pmz) for n,m in pmz_f_temp.iterrows()]
            entropy_denoised_temp = [so.entropy_identity(msms_denoised, m[spectrum_col], pmz) for n,m in pmz_f_temp.iterrows()]
        else:
            pmz_f_temp = pmz_candidate[pmz_candidate['processed_formulas'].isna()]
            entropy_raw_temp = [so.entropy_identity(msms, m[spectrum_col], pmz) for n,m in pmz_f_temp.iterrows()]
            msms_denoised_temp = [raw_msms]*len(pmz_f_temp)
            entropy_denoised_temp = entropy_raw_temp.copy()
        entropy_raw.extend(entropy_raw_temp)
        msms_denoised_all.extend(msms_denoised_temp)
        entropy_denoised.extend(entropy_denoised_temp)
        pmz_candidate_curated = pd.concat([pmz_candidate_curated, pmz_f_temp], ignore_index=True)
    pmz_candidate_curated['msms_denoised']=msms_denoised_all
    pmz_candidate_curated['entropy_raw']=entropy_raw
    pmz_candidate_curated['entropy_denoised']=entropy_denoised
    pmz_candidate_curated.sort_values(by = 'entropy_denoised', ascending = False, inplace = True)
    return(pmz_candidate_curated)
def denoising_search_astral(row, mega_mona):
    raw_msms = (row['spectrum'])
    pmz_candidate = quick_search_sorted(mega_mona, 'PrecursorMZ', row['PrecursorMZ']-0.01, row['PrecursorMZ']+0.01)
    charge = [break_adduct(x)[2] for x in pmz_candidate['Precursor_type']]
    processed_formulas = [prep_formula(x['SMILES'], x['Precursor_type']) for i, x in pmz_candidate.iterrows()]
    pmz_candidate['processed_formulas']=processed_formulas
    pmz_candidate['charge']=charge
    entropy_denoised = []
    entropy_raw = []
    msms_denoised_all = []
    pmz_candidate_curated = pd.DataFrame()
    for f in pmz_candidate['processed_formulas'].unique():

        if f ==f:
            pmz_f_temp = string_search(pmz_candidate, 'processed_formulas', f)
            element_dict, element_count, element_mass, all_possible_mass, all_possible_candidate_formula = drf.get_all_subsets(f)

            msms_denoised = denoising_with_subset(raw_msms, row['PrecursorMZ'], element_dict, all_possible_mass, all_possible_candidate_formula)
            msms_denoised_temp = [msms_denoised]*len(pmz_f_temp)
            entropy_raw_temp = [so.entropy_identity(row['spectrum'], m['spectrum'], row['PrecursorMZ']) for n,m in pmz_f_temp.iterrows()]
            entropy_denoised_temp = [so.entropy_identity(msms_denoised, m['spectrum'], row['PrecursorMZ']) for n,m in pmz_f_temp.iterrows()]
        else:
            pmz_f_temp = pmz_candidate[pmz_candidate['processed_formulas'].isna()]
            entropy_raw_temp = [so.entropy_identity(row['spectrum'], m['spectrum'], row['PrecursorMZ']) for n,m in pmz_f_temp.iterrows()]
            msms_denoised_temp = [raw_msms]*len(pmz_f_temp)
            entropy_denoised_temp = entropy_raw_temp.copy()
        entropy_raw.extend(entropy_raw_temp)
        msms_denoised_all.extend(msms_denoised_temp)
        entropy_denoised.extend(entropy_denoised_temp)
        pmz_candidate_curated = pd.concat([pmz_candidate_curated, pmz_f_temp], ignore_index=True)
    pmz_candidate_curated['msms_denoised']=msms_denoised_all
    pmz_candidate_curated['entropy_raw']=entropy_raw
    pmz_candidate_curated['entropy_denoised']=entropy_denoised
    return(pmz_candidate_curated)
def mass_to_formula(mass, formula, mass_tolerance = 0.01, ifsmile = False):
    if ifsmile == True:
        master_formula = CalcMolFormula(Chem.MolFromSmiles(formula))
    else:
        master_formula = formula
    element_dict, element_count, element_mass, all_possible_mass, all_possible_candidate_formula = get_all_subsets(master_formula)
    left_idx, right_idx =all_possible_mass.searchsorted([mass-mass_tolerance, mass+mass_tolerance+1E-9])
    candidate_subformulas = candidates_to_formulas(all_possible_candidate_formula[left_idx:right_idx], element_dict)
    return candidate_subformulas
def entropy_denoising(msms):
    mass, intensity = so.break_spectra(msms)
    order = np.argsort(intensity)
    mass = mass[order]
    intensity = intensity[order]
    mass_confirmed = np.array([])
    intensity_confirmed = np.array([])
    while len(intensity)>0:
        seed_intensity = np.max(intensity)
        # idx_left = np.searchsorted(intensity, seed_intensity*0.99, side= 'left')
        idx_left = np.searchsorted(intensity, seed_intensity*0.999, side= 'left')
        mass_temp = mass[idx_left:]
        intensity_temp = intensity[idx_left:]
        if len(mass_temp)<=3:
            mass_confirmed =  np.concatenate((mass_confirmed, mass_temp))
            intensity_confirmed = np.concatenate((intensity_confirmed,intensity_temp))
        intensity = intensity[0:idx_left]
        mass = mass[0:idx_left]
    if len(mass_confirmed)==0:
        return np.NAN
    order = np.argsort(mass_confirmed)
    mass_confirmed = mass_confirmed[order]
    intensity_confirmed = intensity_confirmed[order]
    return((so.pack_spectra(mass_confirmed, intensity_confirmed)) )
def dnl_denoising(msms):
    mass, intensity = so.break_spectra(msms)
    order= np.argsort(intensity)
    mass = mass[order]
    intensity = intensity[order]
    from sklearn.linear_model import LinearRegression
    if intensity[1]/2>=intensity[0]*1.5:
        signal_idx = 1

    else:
        k = 2
        if len(mass)==2:
            return(so.pack_spectra([],[]))
        for k in range(2, len(mass)):
            I = intensity[0:k]
            i = np.arange(1,k+1)
            model = LinearRegression().fit(i.reshape((-1,1)), I)
            i_predicted = model.predict(np.array([k+1]).reshape(-1,1))
            if intensity[k]/ i_predicted >2:

                break
        signal_idx = k
    mass_signal = mass[signal_idx:]
    intensity_signal = intensity[signal_idx:]

    return(so.sort_spectrum(so.pack_spectra(mass_signal, intensity_signal)))

def ms_reduce(msms, reduce_factor = 90):
    mass, intensity = so.break_spectra(msms)
    n_chose_peak = np.int32(np.ceil(len(mass)*(1-reduce_factor/100)))
    order = np.argsort(intensity)
    mass = mass[order]
    intensity = intensity[order]
    mass_taken = np.array([])
    intensity_taken = np.array([])
    for i in range(0,11):
        idx_left = np.searchsorted(intensity, np.max(intensity)*(11-i-1)/(11), side = 'left')
        idx_right = np.searchsorted(intensity, np.max(intensity)*(11-i)/(11), side = 'right')

        factor = (n_chose_peak-len(mass_taken))/(idx_right-idx_left)
        if factor>1:
            factor = 1
        sampled_n = np.int32(np.floor(factor*(idx_right-idx_left)))
        sampled_mass = np.random.choice(mass[idx_left:idx_right], size=sampled_n, replace=False)
        sampled_intensity = np.random.choice(intensity[idx_left:idx_right], size=sampled_n, replace=False)
        mass_taken = np.concatenate([mass_taken, sampled_mass])
        intensity_taken= np.concatenate([intensity_taken, sampled_intensity])
        if factor<1:
            break

    return so.sort_spectrum(so.pack_spectra(mass_taken, intensity_taken))
def threshold_denoising(msms, threshold = 1):
    mass, intensity = so.break_spectra(msms)
    intensity_percent = intensity/np.max(intensity)
    to_keep = intensity_percent>(threshold/100)
    mass = mass[to_keep]
    intensity = intensity[to_keep]
    return(so.pack_spectra(mass, intensity))
def spectral_denoising(msms, smiles, adduct,max_allowed_deviation = 0.005):
    # check for empty msms passed, rarely happen
    if isinstance(msms, float):
        return np.NAN
    msms_d1 = entropy_denoising(msms)
    if isinstance(msms_d1, float):
        return np.NAN
    # msms_d1 = entropy_denoising_alt(msms, threshold = 0.001)
    msms_d1_2 = formula_denoising(msms_d1, smiles, adduct,max_allowed_deviation)
    return (msms_d1_2)





def denoising_with_subset(msms, pmz, element_dict, all_possible_mass, all_possible_candidate_formula, mass_tolerance = 0.005):
    msms= entropy_denoising(msms)
    mass, intensity = so.break_spectra(msms)
    parent_ion,mass_tolerance = find_actual_parent_pmz(msms , pmz)
    index_start = np.searchsorted(mass, pmz-1.6,side = 'left')
    mass_parent = mass[index_start:]
    intensity_parent = intensity[index_start:]
    mass_frag = mass[0:index_start]
    intensity_frag = intensity[0:index_start]
    mass_denoised = np.zeros(len(mass_frag))
    intensity_denoised = np.zeros(len(mass_frag))
    counter = 0
    for m,i in zip(mass_frag, intensity_frag):
        left_idx, right_idx =all_possible_mass.searchsorted([parent_ion-m-mass_tolerance, parent_ion-m+mass_tolerance+1E-9])
        candidate_subformulas = candidates_to_formulas(all_possible_candidate_formula[left_idx:right_idx], element_dict)
        if check_candidates(candidate_subformulas):
            mass_denoised[counter]=m
            intensity_denoised[counter]=i
            counter = counter+1
    mass_denoised = mass_denoised[:counter]
    intensity_denoised=intensity_denoised[:counter]
    mass_return = np.concatenate([mass_denoised, mass_parent])
    intensity_return = np.concatenate([intensity_denoised, intensity_parent])
    return(so.pack_spectra(mass_return, intensity_return))
def formula_denoising(msms, smiles, adduct,max_allowed_deviation = 0.005):

    # if type(msms) is list:
    #     msms = so.convert_nist_to_string(msms)
    mol = Chem.MolFromSmiles(smiles)
    # max_c_length = find_longest_element_chain(mol, 'C')
    formula = CalcMolFormula(mol)
    msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    master_formula = prep_formula(smiles, adduct)
    if master_formula!= master_formula:
        return(np.NAN)
    therotical_pmz = calculate_precursormz(formula, adduct)

    index_start = np.searchsorted(mass, therotical_pmz-1.6,side = 'left')
    mass_parent = mass[index_start:]
    intensity_parent = intensity[index_start:]

    if len(mass_parent)==0:
        parent_ion=therotical_pmz
        actual_deviation = 0
    else:
        parent_ion, actual_deviation = find_actual_parent_pmz(msms , therotical_pmz)
    if actual_deviation > max_allowed_deviation:
        mass_tolerance=actual_deviation*1.5
    else:
        mass_tolerance=max_allowed_deviation
    if has_benzene(mol):

        parent_ion = parent_ion+Formula('N2O').isotope.mass
    mass_frag = mass[0:index_start]
    intensity_frag = intensity[0:index_start]
    mass_denoised = np.zeros(len(mass_frag))
    intensity_denoised = np.zeros(len(mass_frag))
    counter = 0
    element_dict, element_count, element_mass, all_possible_mass, all_possible_candidate_formula = get_all_subsets(master_formula)
    for m,i in zip(mass_frag, intensity_frag):
        left_idx, right_idx =all_possible_mass.searchsorted([parent_ion-m-mass_tolerance, parent_ion-m+mass_tolerance+1E-9])
        candidate_subformulas = candidates_to_formulas(all_possible_candidate_formula[left_idx:right_idx], element_dict)
        if check_candidates(candidate_subformulas) == True:
            mass_denoised[counter]=m
            intensity_denoised[counter]=i
            counter = counter+1
    mass_denoised = mass_denoised[:counter]
    intensity_denoised=intensity_denoised[:counter]
    # ei = np.sum(intensity_denoised)/np.sum(intensity_frag)
    mass_return = np.concatenate([mass_denoised, mass_parent])
    intensity_return = np.concatenate([intensity_denoised, intensity_parent])

    return(so.pack_spectra(mass_return, intensity_return))
def get_ei(msms_denoised, msms_raw, pmz):
    if isinstance(msms_denoised, float) and isinstance(msms_raw, float) == False:
        return 0
    msms_denoised = so.truncate_spectrum(msms_denoised, pmz-1.6)
    msms_raw = so.truncate_spectrum(msms_raw, pmz -1.6)
    if isinstance(msms_raw, float) or isinstance(msms_denoised, float):
        return 0

    mass_d, intensity_d = so.break_spectra(msms_denoised)
    mass_r, intensity_r = so.break_spectra(msms_raw)
    return(np.sum(intensity_d)/np.sum(intensity_r))
# def spectra_denoising_frag(msms, smiles, adduct, mass_tolerance = 0.01):
#
#     # max_c_length = find_longest_element_chain(mol, 'C')
#     mol = Chem.MolFromSmiles(smiles)
#     msms = so.sort_spectrum(msms)
#     mass, intensity = so.break_spectra(msms)
#     master_formula = prep_formula_frag(smiles, adduct)
#     therotical_pmz = calculate_precursormz(smiles, adduct, if_smiles=True)
#     index_start = np.searchsorted(mass, therotical_pmz-1.6,side = 'left')
#     mass_parent = mass[index_start:]
#     intensity_parent = intensity[index_start:]
#     mass_frag = mass[0:index_start]
#     intensity_frag = intensity[0:index_start]
#     mass_denoised = np.zeros(len(mass_frag))
#     intensity_denoised = np.zeros(len(mass_frag))
#     counter = 0
#     if adduct[-1]=='+':
#         adjustment = +0.0005500000000040473
#     elif adduct[-1]=='-':
#         adjustment = -0.0005500000000040473
#     else:
#         return(np.NAN, np.NAN)
#     element_dict, element_count, element_mass, all_possible_mass, all_possible_candidate_formula = get_all_subsets(master_formula)
#     for m,i in zip(mass_frag, intensity_frag):
#         left_idx, right_idx =all_possible_mass.searchsorted([m-mass_tolerance+adjustment, m+adjustment+mass_tolerance+1E-9])
#         candidate_subformulas = candidates_to_formulas(all_possible_candidate_formula[left_idx:right_idx], element_dict)
#         if check_candidates_frag(candidate_subformulas) == True:
#             mass_denoised[counter]=mass_denoised[counter]+m
#             intensity_denoised[counter]=intensity_denoised[counter]+i
#             counter = counter+1
#     mass_denoised = mass_denoised[:counter]
#     intensity_denoised=intensity_denoised[:counter]
#     ei = np.sum(intensity_denoised)/np.sum(intensity_frag)
#     mass_return = np.concatenate([mass_denoised, mass_parent])
#     intensity_return = np.concatenate([intensity_denoised, intensity_parent])
#
#     return(so.pack_spectra(mass_return, intensity_return), ei)
def is_denoisable(msms, ref_msms, n = 3):
    msms = so.sort_spectrum(msms)
    mass, intensity = so.break_spectra(msms)
    mass_ref, intensity_ref = so.break_spectra(ref_msms)
    if len(mass_ref)<n or len(mass)<n:
        return False
    order = np.argsort(intensity_ref)
    intensity_ref = intensity_ref[order]
    mass_ref = mass_ref[order]
    ref_key_mass = mass_ref[-n:]
    ref_key_intensity = intensity_ref[-n:]
    ref_key_mass = np.flip(ref_key_mass)
    ref_key_intensity = np.flip(ref_key_intensity)
    key_mass = np.zeros(n)
    key_intensity = np.zeros(n)
    ratios = np.zeros(n-1)
    ref_ratios = np.zeros(n-1)
    for i in range(0, n):
        idx_left, idx_right = mass.searchsorted([ref_key_mass[i]-0.01, ref_key_mass[i]+0.01+1E-9])
        if idx_right-idx_left==0:
            return( False ) #one of the key mass is not found
        anchor = np.argmax(intensity[idx_left:idx_right])
        key_mass[i]=mass[idx_left:idx_right][anchor]
        key_intensity[i]=intensity[idx_left:idx_right][anchor]
    for i in range(0, n-1):
        ratios[i]=key_intensity[i]/key_intensity[i+1]
        ref_ratios[i]=ref_key_intensity[i]/ref_key_intensity[i+1]
    cross_ratio = ref_ratios/ratios
    if np.max(np.diff(key_intensity))<0:
        return True
    else:
        return False
# def check_candidates_frag(candidates):
#     if len(candidates)==0:
#         return False
#     for c in candidates:
#         if check_candidate_frag(c)==True:
#             return True
#     return False
# def check_candidate_frag(subformula):
#     if 'H' in subformula:
#         candidate_variations = [Formula(subformula).__sub__(Formula('H')).formula, subformula,Formula(subformula).__add__(Formula('H')).formula]
#     else:
#         candidate_variations = [ subformula,Formula(subformula).__add__(Formula('H')).formula]
#     for c in candidate_variations:
#         # if check_senior(c) and check_ratio(c)==True:
#         if check_ratio(c)==True and check_senior(c):
#             return True
#     return False
def find_actual_parent_pmz(msms,therotcial_pmz):
    mass_parent, intensity_parent = so.break_spectra(msms)
    if len(mass_parent)==0:
        return(therotcial_pmz)
    left_idx, right_idx = mass_parent.searchsorted([therotcial_pmz-0.01, therotcial_pmz+0.01])
    if right_idx-left_idx==0:
        return therotcial_pmz, 0.01
    max_idx = np.argmax(intensity_parent[left_idx:right_idx])
    return(mass_parent[left_idx:right_idx][max_idx], abs(mass_parent[left_idx:right_idx][max_idx]-therotcial_pmz)*1.2)

def check_single_candidate(subformula):
    return check_ratio(subformula)
    # if 'H' in subformula:
    #     candidate_variations = [Formula(subformula).__sub__(Formula('H')).formula, subformula,Formula(subformula).__add__(Formula('H')).formula]
    # else:
    #     candidate_variations = [ subformula,Formula(subformula).__add__(Formula('H')).formula]
    # for c in candidate_variations:
    #     # if check_senior(c) and check_ratio(c)==True:
    #     if check_ratio(c)==True:
    #         return True
    # return False

def check_candidates(candidate_subformulas):
    if len(candidate_subformulas)==0:
        return False
    for c in candidate_subformulas:
        if check_single_candidate(c)==True:
            return True
    return False
def get_all_subsets(master_formula):
    master_formula = chemparse.parse_formula(master_formula)
    element_dict=np.fromiter(master_formula.keys(), dtype='U2')
    element_count=np.fromiter(master_formula.values(), dtype=np.int16)
    element_mass = np.array([Formula(x).isotope.mass for x in element_dict])
    formula_range = [range(x + 1) for x in element_count]
    all_possible_candidate_formula = np.array(
        list(itertools.product(*formula_range)), np.int16)
    all_possible_mass = np.sum(
        element_mass * all_possible_candidate_formula, axis=1)
    order = np.argsort(all_possible_mass)
    all_possible_mass = all_possible_mass[order]
    all_possible_candidate_formula = all_possible_candidate_formula[order]
    return (element_dict, element_count, element_mass, all_possible_mass, all_possible_candidate_formula)
from toolsets.constants import valence_dict
def check_senior(formula):
    max_valence = 0
    sum_valence = 0
    element_count = 0
    parsed_formula = chemparse.parse_formula(formula)
    for k in parsed_formula.keys():
        sum_valence= sum_valence+parsed_formula[k]*valence_dict[k]
        element_count = element_count+parsed_formula[k]
        if valence_dict[k]>max_valence:
            max_valence=valence_dict[k]
    if sum_valence%2 != 0:
        return False
    if sum_valence<2*max_valence:
        return False
    if sum_valence<2*(element_count-1):
        return False
    return True

def check_ratio(formula):
    if len(formula)==0:
        return False
    parsed_formula = chemparse.parse_formula(formula)
    element_dict=np.fromiter(parsed_formula.keys(), dtype='U2')
    element_count=np.fromiter(parsed_formula.values(), dtype=np.int16)
    accurate_mass = Formula(formula).isotope.mass
    if accurate_mass >0 and accurate_mass<500:
        if 'H' in element_dict and parsed_formula['H']>72:
            return False
        if 'C' in element_dict and parsed_formula['C']>39:
            return False
    elif accurate_mass>500 and accurate_mass<1000:
        if 'H' in element_dict and parsed_formula['H']>126:
            return False
        if 'C' in element_dict and parsed_formula['C']>78:
            return False

    if 'N' in element_dict and parsed_formula['N']>20:
        return False

    if 'O' in element_dict and parsed_formula['O']>27:
        return False
    if len(element_dict)==1 and element_dict[0]=='C':
        #check for pure carbon loss
        return(False)
    if len(element_dict)==1 and element_dict[0]=='N' and element_count[0]!=2:
        #check for pure nitrogen loss (while not N2)
        return (False)
    if 'C' in element_dict and 'H' in element_dict:
        if parsed_formula['H']/parsed_formula['C']>6 or parsed_formula['H']/parsed_formula['C']<0.1:
            # 7 golden rules: HC check
            return False
    if 'C' in element_dict and 'F' in element_dict:
        if parsed_formula['F']/parsed_formula['C']>6:
            # 7 golden rules: CF check
            return False

    if 'C' in element_dict and 'Cl' in element_dict:
        if parsed_formula['Cl']/parsed_formula['C']>2:
            # 7 golden rules: CCl check
            return False

    if 'C' in element_dict and 'Br' in element_dict:
        if parsed_formula['Br']/parsed_formula['C']>2:
            # 7 golden rules: CF check
            return False

    if 'C' in element_dict and 'N' in element_dict:
        if parsed_formula['N']/parsed_formula['C']>4:
            # 7 golden rules: CF check
            return False

    if 'C' in element_dict and 'O' in element_dict:
        if parsed_formula['O']/parsed_formula['C']>3:
            # 7 golden rules: CF check
            return False

    if 'C' in element_dict and 'P' in element_dict:
        if parsed_formula['P']/parsed_formula['C']>2:
            # 7 golden rules: CF check
            return False

    if 'C' in element_dict and 'S' in element_dict:
        if parsed_formula['S']/parsed_formula['C']>3:
            # 7 golden rules: CF check
            return False

    if 'C' in element_dict and 'Si' in element_dict:
        if parsed_formula['Si']/parsed_formula['C']>1:
            # 7 golden rules: CF check
            return False
    return(True)
# def prep_formula_frag(smiles, adduct):
#
#     mol = Chem.MolFromSmiles(smiles)
#     extra_atoms = has_benzene(mol)
#     # max_c_length = find_longest_element_chain(mol, 'C')
#
#     formula = CalcMolFormula(mol)
#
#     if adduct in ['[M]+', '[Cat]+', '[M]-']:
#         if formula[-1] in ['+', '-'] and formula[-1]==adduct[-1]:
#             formula = formula[0:-1]
#             master_formula = Formula(formula)
#             if extra_atoms ==True:
#                 master_formula = master_formula.__add__(Formula('N2O'))
#             return(master_formula)
#         else:
#             return(np.NAN)
#
#     if 'Hac' in adduct:
#         adduct = adduct.replace("Hac", "C2H4O2")
#     if 'FA' in adduct:
#         adduct = adduct.replace("FA", "CH2O2")
#     if 'DMSO' in adduct:
#         adduct = adduct.replace("DMSO", "C2H6OS")
#     if adduct[0]=='[':
#         adduct_part = re.split(r'\[|]', adduct)[1]
#         ind = re.split(r'(\+|-)', adduct_part)
#     else:
#         adduct_part = adduct
#         ind = re.split(r'(\+|-)', adduct_part)
#     master_formula = Formula(formula)
#     coef,a  = parse_adduct(ind[0])
#     for i in range(1, coef):
#         master_formula = master_formula.__add__(master_formula)
#     for i in range(1, len(ind)):
#         if (ind[i] not in ['+', '-']) and len(ind[i])>0:
#             coef, a = parse_adduct(ind[i])
#             if ind[i-1]=='+':
#                 for j in range(0, coef):
#                     master_formula = master_formula.__add__(Formula(a))
#             elif ind[i-1]=='-':
#                 for j in range(0, coef):
#                     master_formula = master_formula.__sub__(Formula(a))
#     if extra_atoms == True:
#         master_formula = master_formula.__add__(Formula('N2O'))
#     return(master_formula.formula)
def prep_formula(smiles, adduct, show_original = False):
    if smiles != smiles or adduct != adduct or 'i' in adduct:
        return np.NAN
    mol = Chem.MolFromSmiles(smiles)
    if adduct not in ['[M]+', '[Cat]+', 'M+'] and adduct[-1]=='+' and Chem.GetFormalCharge(mol)<0:
        un = rdMolStandardize.Uncharger()
        mol = un.uncharge((mol))
        if Chem.GetFormalCharge(mol)!=0:
            return(np.NAN)


    extra_atoms = has_benzene(mol)
    if np.abs(Chem.GetFormalCharge(mol))>1:
        return(np.NAN)
        # max_c_length = find_longest_element_chain(mol, 'C')

    formula = CalcMolFormula(mol)
    if show_original == True:
        print(formula)
    if adduct in ['[M]+', '[Cat]+', '[M]-']:
        if formula[-1] in ['+', '-'] and formula[-1]==adduct[-1]:
            formula = formula[0:-1]
            master_formula = Formula(formula)
            if extra_atoms ==True:
                master_formula = master_formula.__add__(Formula('N2O'))
            return(master_formula.formula)
        else:
            return(np.NAN)

    if 'Hac' in adduct:
        adduct = adduct.replace("Hac", "C2H4O2")
    if 'FA' in adduct:
        adduct = adduct.replace("FA", "CH2O2")
    if 'DMSO' in adduct:
        adduct = adduct.replace("DMSO", "C2H6OS")
    adduct_part, ind, charge = break_adduct(adduct)
    if abs(charge) != 1:
        return np.NAN
    master_formula = Formula(formula)
    coef,a  = parse_ind(ind[0])
    for i in range(1, coef):
        master_formula = master_formula.__add__(master_formula)
    for i in range(1, len(ind)):
        if (ind[i] not in ['+', '-', 'H']) and len(ind[i])>0:
            coef, a = parse_ind(ind[i])
            if ind[i-1]=='+':
                for j in range(0, coef):
                    master_formula = master_formula.__add__(Formula(a))
            elif ind[i-1]=='-':
                for j in range(0, coef):
                    master_formula = master_formula.__sub__(Formula(a))
    if extra_atoms == True:
        master_formula = master_formula.__add__(Formula('N2O'))

    return(master_formula.formula)

def has_benzene(molecule):
    # molecule = Chem.MolFromSmiles(smiles)
    benzene = Chem.MolFromSmiles('c1ccccc1')  # Aromatic benzene SMILES notation

    # Check if benzene is a substructure of the given molecule
    return molecule.HasSubstructMatch(benzene)
def candidates_to_formulas(candidates, element_dict):
    formulas = []
    for c in candidates:
        formulas.append(dict_to_formula(c, element_dict))
    return formulas
def dict_to_formula(candidate, element_dict):
    string = ''
    for i in range(0, len(candidate)):
        if candidate[i]>1:
            string += element_dict[i] + str(candidate[i])
        elif candidate[i]==1:
            string += element_dict[i]
    return string
#
#
# def denoise_h(msms, smiles, adduct, parent_ion,reference_db_sorted, mass_error = 0.01):
#     start = time.time()
#     msms = so.sort_spectrum(msms)
#     mass, intensity = so.break_spectra(msms)
#
#     index_start = np.searchsorted(mass, parent_ion-1.6,side = 'left')
#     mass_parent = mass[index_start:]
#     intensity_parent = intensity[index_start:]
#     mass_frag = mass[0:index_start]
#     intensity_frag = intensity[0:index_start]
#     formula = prep_formula(smiles, adduct)
#     mass_d = []
#     intensity_d = []
#     mass_tbd = []
#     intensity_tbd = []
#     master_formula = MolecularFormula()
#     master_formula.from_string(formula)
#     for cur_i in range(0,len(mass_frag)):
#         loss_candidate_white = get_loss_candidate(parent_ion-mass_frag[cur_i], reference_db_sorted, step = mass_error)
#         if check_loss(loss_candidate_white, master_formula)==True:
#             mass_d.append(mass_frag[cur_i])
#             intensity_d.append(intensity_frag[cur_i])
#         else:
#             mass_tbd.append(mass_frag[cur_i])
#             intensity_tbd.append(intensity_frag[cur_i])
#     master_formula = MolecularFormula()
#     master_formula.from_string(formula)
#     precursor_data = master_formula.get_data()
#     formula_range = [range(x + 1) for x in precursor_data]
#     all_possible_candidate_formula = np.array(
#         list(itertools.product(*formula_range)), numpy_formula_format)
#     all_possible_mass = np.sum(
#         atom_mass_array * all_possible_candidate_formula, axis=1)
#     order = np.argsort(all_possible_mass)
#     all_possible_mass = all_possible_mass[order]
#     all_possible_candidate_formula = all_possible_candidate_formula[order]
#
#     for i in range(len(mass_tbd)):
#         loss_mass = parent_ion - mass_tbd[i]
#         if check_frag_bl(loss_mass, all_possible_mass,all_possible_candidate_formula, mass_error):
#             mass_d.append(mass_tbd[i])
#             intensity_d.append(intensity_tbd[i])
#     # ei = so.calculate_explained_intensity(so.pack_spectra(mass_d, intensity_d), so.pack_spectra(mass_frag, intensity_frag))
#     if len(mass_d)==0:
#         return(np.NAN, 0)
#     ei = np.sum(intensity_d)/np.sum(intensity_frag)*100
#     mass_d.extend(mass_parent)
#     intensity_d.extend(intensity_parent)
#     return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)), ei)
# def check_frag_bl(loss_mass,all_possible_mass, all_possible_candidate_formula, mass_error):
#     l, r = all_possible_mass.searchsorted([loss_mass-mass_error, loss_mass+mass_error+1E-9])
#     if r-l <1:
#         return False
#
#     for f in all_possible_candidate_formula[l:r]:
#         if evaluate_nl_blacklist(MolecularFormula(f).__str__())==True:
#             return(True)
#     return(False)
# def check_loss(loss_formula_candidate, master_formula):
#     # print(formula)
#     for c in loss_formula_candidate:
#         # candidate_formula = mtf.MolecularFormula()
#         # try:
#         #     candidate_formula.from_string(c)
#         # except:
#         #     continue
#         if check_candidate(c, master_formula)==True:
#             return True
#
#     return(False)
# def check_candidate(c, master_formula):
#     candidate_formula = MolecularFormula()
#     try:
#         candidate_formula.from_string(c)
#     except:
#         return(False)
#     for i in range(0, len(master_formula.get_data())):
#         if candidate_formula.get_data()[i]>master_formula.get_data()[i]:
#             return(False)
#     return(True)
# def evaluate_nl_blacklist(nl):
#     if len(nl)==0:
#         return False
#
#     mol = MolecularFormula()
#     mol.from_string(nl)
#     if np.max(np.append(mol.get_data()[0:1], mol.get_data()[2:]))==0 or np.max(np.append(mol.get_data()[0:2], mol.get_data()[3:]))==0:
#
#         return False
#     n_H = mol.get_data()[atom_list.index('H')]
#     n_C = mol.get_data()[atom_list.index('C')]
#     n_N = mol.get_data()[atom_list.index('N')]
#     n_O = mol.get_data()[atom_list.index('O')]
#     n_F = mol.get_data()[atom_list.index('F')]
#     n_Na = mol.get_data()[atom_list.index('Na')]
#     n_P = mol.get_data()[atom_list.index('P')]
#     n_S = mol.get_data()[atom_list.index('S')]
#     n_Cl = mol.get_data()[atom_list.index('Cl')]
#     n_K = mol.get_data()[atom_list.index('K')]
#     n_Br = mol.get_data()[atom_list.index('Br')]
#     n_I = mol.get_data()[atom_list.index('I')]
#     if n_C != 0 or n_O != 0 or n_N != 0 or n_P != 0 or n_S != 0:
#         if n_H+n_Br+n_I+n_F+n_Cl>n_C*4+n_O*4+n_N*5+n_P*5+n_S*4:
#
#             return False
#     if n_Na >1:
#         return(False)
#     if n_N/n_C >4 and n_C>0:
#         return(False)
#     if n_O/n_C>3 and n_C>0:
#         return(False)
#     if n_H/n_C>6 and n_C>0:
#         return(False)
#     if n_P/n_C>2 and n_C>0:
#         return False
#     if n_K >1:
#         return(False)
#     if n_Na >0 and n_H >0:
#         if np.max(np.append(mol.get_data()[1:5], mol.get_data()[6:]))==0:
#             return False
#     if n_K >0 and n_H >0:
#         if np.max(np.append(mol.get_data()[1:9], mol.get_data()[10:]))==0:
#             return False
#     return True
# def get_loss_candidate(mass, reference_db, step):
#     mass_neg = mass- atom_mass['H']
#     mass_pos = mass+atom_mass['H']
#     losses = []
#     search_array = np.array(reference_db['mass'])
#     for m in [mass_neg, mass, mass_pos]:
#         # print(m)
#         l, r = search_array.searchsorted([m-step, m+step+1E-9])
#         if r-l >=1:
#             losses.extend(reference_db.iloc[l:r]['formstr'].tolist())
#         # losses = losses.append(quick_search_sorted(reference_db, 'mass', m, step))
#     return (losses)


# def denoise_h_formula(msms, formula, parent_ion, adduct, reference_db_sorted, mass_error = 0.005):
#     start = time.time()
#     msms = so.sort_spectrum(msms)
#     mass, intensity = so.break_spectra(msms)
#     # parent_ion = calculate_precursormz(smiles, adduct)
#     # precursor = get_precursor(msms, parent_ion)
#     # if precursor == precursor:
#     #     parent_ion = precursor
#     # state = 'inferred'
#     # mass_error = abs(parent_ion-precursor)*3
#     # else:
#     # state = 'default'
#     # mass_error = 0.005
#     index_start = np.searchsorted(mass, parent_ion-1.5,side = 'left')
#     mass_parent = mass[index_start:]
#     intensity_parent = intensity[index_start:]
#     mass_frag = mass[0:index_start]
#     intensity_frag = intensity[0:index_start]
#     # formula = prep_formula(smiles, adduct)
#     if adduct=='[M+Na]+':
#         formula = formula+'Na'
#     elif adduct=='[M+NH4]+':
#         formula = formula+'NH3'
#     elif adduct=='[M+Cl]-':
#         formula = formula+'Cl'
#     elif adduct=='[M+Hac-H]-':
#         formula = formula+'C2H4O2'
#     elif adduct =='[M+K]+':
#         formula = formula+'K'
#     elif adduct == '[M+FA-H]':
#         formula = formula +'CH2O2'
#     else:
#         formula = formula
#     formula = (Formula(formula).formula)
#     mass_d = []
#     intensity_d = []
#     mass_tbd = []
#     intensity_tbd = []
#     for cur_i in range(len(mass_frag)):
#         loss_candidate_white = get_loss_candidate(parent_ion-mass_frag[cur_i], reference_db_sorted, step = mass_error)
#         if (check_loss(loss_candidate_white['formstr'].values.tolist(), formula))==True:
#             mass_d.append(mass_frag[cur_i])
#             intensity_d.append(intensity_frag[cur_i])
#         else:
#             mass_tbd.append(mass_frag[cur_i])
#             intensity_tbd.append(intensity_frag[cur_i])
#     check_1 = time.time()
#     # all_possible_formulas = get_all_formulas(formula)
#     # all_allowed_formulas, all_allowed_masses = get_all_allowed_formula_mass(all_possible_formulas)
#     for i in range(len(mass_tbd)):
#         loss_mass = parent_ion - mass_tbd[i]
#         if check_frag_bl(loss_mass, formula, mass_error):
#             mass_d.append(mass_tbd[i])
#             intensity_d.append(intensity_tbd[i])
#     ei = so.calculate_explained_intensity(so.pack_spectra(mass_d, intensity_d), so.pack_spectra(mass_frag, intensity_frag))
#     mass_d.extend(mass_parent)
#     intensity_d.extend(intensity_parent)
#     return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)), ei)
#     # return(so.pack_spectra(mass_d, intensity_d), state, mass_error)
# def denoise_w(msms, smiles, adduct, reference_db_sorted, mass_error = 0.005):
#     start = time.time()
#     msms = so.sort_spectrum(msms)
#     mass, intensity = so.break_spectra(msms)
#     parent_ion = calculate_precursormz(smiles, adduct)
#     # precursor = get_precursor(msms, parent_ion)
#     # if precursor == precursor:
#     #     parent_ion = precursor
#     # state = 'inferred'
#     # mass_error = abs(parent_ion-precursor)*3
#     # else:
#     # state = 'default'
#     # mass_error = 0.005
#     index_start = np.searchsorted(mass, parent_ion-1.5,side = 'left')
#     mass_parent = mass[index_start:]
#     intensity_parent = intensity[index_start:]
#     mass_frag = mass[0:index_start]
#     intensity_frag = intensity[0:index_start]
#     formula = prep_formula(smiles, adduct)
#     mass_d = []
#     intensity_d = []
#     mass_tbd = []
#     intensity_tbd = []
#     for cur_i in range(len(mass_frag)):
#         loss_candidate_white = get_loss_candidate(parent_ion-mass_frag[cur_i], reference_db_sorted, step = mass_error)
#         if (check_loss(loss_candidate_white['formstr'].values.tolist(), formula))==True:
#             mass_d.append(mass_frag[cur_i])
#             intensity_d.append(intensity_frag[cur_i])
#         else:
#             mass_tbd.append(mass_frag[cur_i])
#             intensity_tbd.append(intensity_frag[cur_i])
#     check_1 = time.time()
#     # all_possible_formulas = get_all_formulas(formula)
#     # all_allowed_formulas, all_allowed_masses = get_all_allowed_formula_mass(all_possible_formulas)
#     # for i in range(len(mass_tbd)):
#     #     loss_mass = parent_ion - mass_tbd[i]
#     #     if check_frag_bl(loss_mass, formula, mass_error):
#     #         mass_d.append(mass_tbd[i])
#     #         intensity_d.append(intensity_tbd[i])
#     ei = so.calculate_explained_intensity(so.pack_spectra(mass_d, intensity_d), so.pack_spectra(mass_frag, intensity_frag))
#     mass_d.extend(mass_parent)
#     intensity_d.extend(intensity_parent)
#     return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)), ei)
# def check_frag_bl(frag_mass, formula, mass_error = 0.005):
#     nl_candidates = mtf.nl_to_formula(frag_mass, mass_error, formula)
#     if nl_candidates != [] and nl_candidates!=['']:
#         for nl in nl_candidates:
#             if evaluate_nl_blacklist(nl)==True:
#                 return(True)
#         else:
#             return(False)
#     else:
#         return(False)
# # if intensity_precursor[0] != -1:
# #     mass_d.extend(mass_precursor)
# #     intensity_d.extend(intensity_precursor)
#
# def find_loss(loss_mass, mass_error, all_allowed_formulas, all_allowed_masses):
#     index_start = np.searchsorted(all_allowed_masses, loss_mass-mass_error,side = 'left')
#     index_end = np.searchsorted(all_allowed_masses, loss_mass+mass_error,side = 'right')
#     if len(all_allowed_masses[index_start:index_end])>0:
#         return(all_allowed_formulas[index_start:index_end])
#     else:
#         return ([])
# def get_all_allowed_formula_mass(all_possible_formulas):
#     all_allowed_formulas = []
#     all_allowed_masses = []
#     check_21 = time.time()
#     for formula in all_possible_formulas:
#         if evaluate_nl_blacklist(formula)==True:
#             all_allowed_formulas.append(formula)
#             all_allowed_masses.append(Formula(formula).isotope.mass)
#     if len(all_allowed_formulas)<1:
#         return([], [])
#     allowed_mass_sorted, allowed_formula_sorted = zip(*sorted(zip(all_allowed_masses, all_allowed_formulas)))
#     return(allowed_formula_sorted,allowed_mass_sorted )
# def get_all_formulas(formula):
#     precursor_formula = MolecularFormula()
#     precursor_formula.from_string(Formula(formula).formula)
#     precursor_data = precursor_formula.get_data()
#     formula_range = [range(x + 1) for x in precursor_data]
#     all_possible_candidate_formula = np.array(
#         list(itertools.product(*formula_range)), numpy_formula_format)
#     all_possible_formulas = translate_matrix_into_formula(all_possible_candidate_formula)
#     return(all_possible_formulas)
# def translate_matrix_into_formula(matrix):
#     formulas = []
#     for i in matrix:
#         if np.max(i)>0:
#             formula_temp = ''
#             for j in range(len(i)):
#                 if i[j]>0:
#                     if i[j]>1:
#                         formula_temp = formula_temp+atom_list[j]+str(i[j])
#                     elif i[j]==1:
#                         formula_temp = formula_temp+atom_list[j]
#             formulas.append(formula_temp)
#     return formulas
# def get_precursor(msms, precursor_mz_tho):
#     msms = so.sort_spectrum(msms)
#     mass, intensity = so.break_spectra(msms)
#     index_start = np.searchsorted(mass, precursor_mz_tho-1.5,side = 'left')
#     mass_parent = mass[index_start:]
#     intensity_parent = intensity[index_start:]
#     if len(mass_parent)>0:
#         precursor_idx = np.argmax(intensity_parent)
#         return(mass_parent[precursor_idx])
#     else:
#         return(np.NAN)
#
# from rdkit.Chem.AllChem import CalcExactMolWt
# from rdkit import Chem
# from toolsets.constants import single_charged_adduct_mass
# def denoise_hybrid_msp(instance, reference_db_sorted, typeofmsms = 'peaks_recalibrated', mass_error =0.02 ):
#     mass, intensity = so.break_spectra(instance[typeofmsms])
#     adduct = instance['Precursor_type']
#     # parention = float(instance['PrecursorMZ'])
    #     parention = CalcExactMolWt(Chem.MolFromSmiles(instance['SMILES']))+single_charged_adduct_mass[adduct]
#
#     formula = prep_formula(instance['SMILES'], adduct)
#     mass_d = []
#     intensity_d = []
#     for cur_i in range(len(mass)):
#         loss_candidate_white = get_loss_candidate(parention-mass[cur_i], reference_db_sorted, step = mass_error)
#         if (check_loss(loss_candidate_white['formstr'].values.tolist(), formula))==True:
#             mass_d.append(mass[cur_i])
#             intensity_d.append(intensity[cur_i])
#         else:
#             loss_candidate_black = mtf.nl_to_formula(parention-mass[cur_i], mass_error, formula)
#             if(check_candidate_black(loss_candidate_black))==True:
#                 mass_d.append(mass[cur_i])
#                 intensity_d.append(intensity[cur_i])
#     return(so.pack_spectra(mass_d, intensity_d))
#
#
#
# def check_candidate_black(loss_candidate_black):
#     if len(loss_candidate_black) > 0 and loss_candidate_black!=['']:
#         for nl in loss_candidate_black:
#             if evaluate_nl_blacklist(nl)==True:
#                 return(True)
#         return(False)
#     else:
#         return(False)

#

# # def check_candidate(loss_temp_dict, formula_dict):
# # #     both loss_temp_d and formula_d are in dataframe form
# #     if set(loss_temp_dict.keys()).issubset(set(formula_dict.keys())):
# #         # is all loss atom subset of formula atom?
# #         for element in loss_temp_dict.keys():
# #             if loss_temp_dict[element] > formula_dict[element]:
# #                 return( False)
# #         return(True)
# #     else:
# #         return(False)
# def denoise_blacklist(instance, typeofmsms = 'peaks_recalibrated', mass_error =0.02):
#     # print("i am in new")
#     import toolsets.spectra_operations as so
#     mass, intensity = so.break_spectra(instance[typeofmsms])
#     parention = float(instance['reference_precursor_mz'])
#     # pep_mass, pep_intensity = get_parentpeak(mass, intensity, parention)
#     adduct = instance['reference_adduct']
#     formula = prep_formula(instance['reference_smiles'], adduct)
#     # mass_frag, intensity_frag, mass_precursor, intensity_precursor = remove_precursor(mass, intensity, parention)
#
#     # if instance['reference_adduct']=='[M-H2O-H]-' or instance['reference_adduct']=='[M+H2O+H]+':
#     #     formula = formula+'H2O'
#     mass_d = []
#     intensity_d = []
#
#     # print(threshold)
#     for i in range(0,len(mass)):
#         # print("i am evaluating", mass[i])
#         nl_candidates = mtf.nl_to_formula(parention - mass[i], mass_error, formula)
#         if nl_candidates != [] and nl_candidates!=['']:
#             for nl in nl_candidates:
#                 if evaluate_nl_blacklist(nl)==True:
#                     mass_d.append(mass[i])
#                     intensity_d.append(intensity[i])
#                     break
#     # if intensity_precursor[0] != -1:
#     #     mass_d.extend(mass_precursor)
#     #     intensity_d.extend(intensity_precursor)
#     return(so.pack_spectra(mass_d, intensity_d))
# def denoise_single(frag_ion, parention, smiles, mass_error = 0.02, adduct = '[M+H]+',ifppm = False, ifnl = True):
#     import toolsets.spectra_operations as so
#     parention = float(parention)
#     threshold = so.set_tolerance(mass_error = mass_error, ifppm=ifppm, precursormz=parention)
#     # formula= str(formula)
#     formula = prep_formula(smiles, adduct)
#     if ifnl ==True:
#         nl_candidates = mtf.nl_to_formula(parention - frag_ion, threshold, formula)
#         return(nl_candidates)
#     if ifnl == False:
#         nl_candidates = mtf.nl_to_formula(frag_ion, threshold, formula)
#         return(nl_candidates)
#

# def post_processing(data, raw_column = 'peaks_recalibrated',denoised_column = 'peaks_denoised'):
#     import toolsets.spectra_operations as so
#     ei = []
#     for index, row in data.iterrows():
#         ei.append(so.calculate_explained_intensity(row[raw_column], row[denoised_column]))
#     data['ei']=ei
#     peaks_denoisied_normalized = []
#     normalized_entropy = []
#     spectrum_entropy = []
#     # data['c_id'] = np.arange(len(data))
#     for index, row in data.iterrows():
#         peaks_denoisied_normalized.append(so.normalize_spectrum(row[denoised_column]))
#         normalized_entropy.append(so.normalized_entropy(row[denoised_column], order=4))
#         spectrum_entropy.append(so.spectral_entropy(row[denoised_column]))
#     data['peaks_denoised_normalized']=peaks_denoisied_normalized
#     data['spectrum_entropy']=spectrum_entropy
#     data['normalized_entropy']=normalized_entropy
#     return(data)
#
#
# from toolsets.mass_to_formula import atom_list
# from toolsets.mass_to_formula import MolecularFormula, atom_mass_array, numpy_formula_format,atom_list
# def get_allowed_formula(all_possible_candidate_formula, element):
#
#     idx = atom_list.index(element)
#     all_allowed_candidate_formula = np.copy(all_possible_candidate_formula)
#     pop_idx = []
#     for i in range(len(all_possible_candidate_formula)):
#         if np.max(all_possible_candidate_formula[i][0:idx])==0 and np.max(all_possible_candidate_formula[i][idx+1:])==0:
#             pop_idx.append(i)
#     all_allowed_candidate_formula = np.delete(all_allowed_candidate_formula, pop_idx, axis=0)
#     # all_allowed_candidate_formula.append(all_possible_candidate_formula[i])
#     return(all_allowed_candidate_formula)
# def get_all_allowed_formula(all_possible_candidate_formula, elements):
#     all_allowed_candidate_formula = np.copy(all_possible_candidate_formula)
#     for element in elements:
#         all_allowed_candidate_formula = get_allowed_formula(all_allowed_candidate_formula, element = element)
#     return(all_allowed_candidate_formula)
# def math_check(smiles, adduct, return_mass_list = False):
#     precursor_mz = calculate_precursormz(smiles, adduct)
#     formula = prep_formula(smiles, adduct)
#     precursor_formula = MolecularFormula()
#
#     precursor_formula.from_string(Formula(formula).formula)
#     precursor_data = precursor_formula.get_data()
#     formula_range = [range(x + 1) for x in precursor_data]
#     forbidden_elements = ['C','N']
#     if 'Na' in formula:
#         forbidden_elements.append('Na')
#     elif 'K' in formula:
#         forbidden_elements.append('K')
#     elif 'Cl' in formula:
#         forbidden_elements.append('Cl')
#     all_possible_candidate_formula = np.array(
#         list(itertools.product(*formula_range)), numpy_formula_format)
#     all_allowed_candidate_formula = get_all_allowed_formula(all_possible_candidate_formula, elements=forbidden_elements)
#     all_possible_mass = np.sum(
#         atom_mass_array * all_allowed_candidate_formula, axis=1)
#     all_possible_mass.sort()
#     diff = np.diff(all_possible_mass)
#     allowed_length = len(diff[diff>0.005])*0.005*2
#     # if return_mass_list == True:
#     #     return(all_possible_mass)
#     # peak_list = make_peak_list(all_possible_mass)
#     # peak_list_checked = check(peak_list)
#     # allowed_region = 0
#     # for i in range(0, len(peak_list_checked)):
#     #     allowed_region = allowed_region+(peak_list_checked[i][2]-peak_list_checked[i][0])
#     return(allowed_length/precursor_mz*100)
# def make_peak_list(all_mass_list):
#     peak_list = []
#     all_mass_list.sort()
#     for mass in all_mass_list:
#         peak_list.append([mass-0.005, mass, mass+0.005])
#     return peak_list
# def check_overlap(peak_list):
#     # peak_list.sort()
#     for p in range(0, len(peak_list)-1):
#         if peak_list[p][2]>peak_list[p+1][0]:
#             return(p)
#     return(-1)
# def check(peak_list):
#
#     while check_overlap(peak_list)!= -1:
#         idx = check_overlap(peak_list)
#         peak_list[idx][2]=peak_list[idx+1][2]
#         peak_list = np.delete(peak_list, idx+1, axis=0)
#
#     return(peak_list)
#     # peak_list[p][2]=peak_list[p+1][2]
#     # data_good=num_search(data, 'ei', ei_threshold, direction=">", inclusion=True)
#     #     # data=num_search(data, 'ei', ei_threshold, direction="<", inclusion=False)
#     # data_good =num_search(data_good, 'spectrum_entropy', 0.5, direction='>', inclusion=True)
#     # data_good.reset_index(drop = True, inplace= True)
#     # if high_quality ==True:
#     #     return(data_good)
#     # else:
#     #     # data['c_id'] = np.arange(len(data_HCD))
#     #     data_bad = data[~data['library_id'].isin(data_good['library_id'])]
#     #     # data_bad =num_search(data_bad, 'spectrum_entropy', 0.5, direction='>', inclusion=True)
#     #     data_bad.reset_index(drop = True, inplace= True)
#     #
#     #     return(data_bad)
# # def prep_formula(formula, adduct = ['[M+H]+']):
# #     formula= str(formula)
# #     if adduct=='[M+Na]+':
# #         formula = formula+'Na'
# #     if adduct=='[M+NH4+]':
# #         formula = formula+'NH3'
# #     if adduct=='[M+Cl]-':
# #         formula = formula+'Cl'
# #     if adduct=='[M+C2H4O2-H]-':
# #         formula = formula+'C2H4O2'
# #     if formula[-1]=='+' or formula[-1]=='-':
# #         formula = formula[:-1]
# #     return(formula)
# # def find_parention(mass, intensity, precursormz, mass_error = 10, ifppm = True):
# #     precursormz = float(precursormz)
# #     tol = so.set_tolerance(mass_error, ifppm, precursormz)
# #     # mass_raw, intensity_raw = so.break_spectra(msms)
# #     # if precursormz <=200:
# #     #     tolerance = 0.01
# #     # if precursormz >200:
# #     #     tolerance = precursormz*tolerance/1E6
# #
# #     # mass_sorted, intensity_sorted = zip(*sorted(zip(mass_raw, intensity_raw)))
# #     lower_bound = precursormz-tol
# #     upper_bound = precursormz+tol
# #     # print(" i am here!")
# #     lower_bound_i = bisect.bisect_left(mass, lower_bound)
# #     # print("i have passed lower")
# #     upper_bound_i = bisect.bisect_right(mass, upper_bound, lo=lower_bound_i)
# #     mass_candi = mass[lower_bound_i:upper_bound_i]
# #     intensity_candi = intensity[lower_bound_i:upper_bound_i]
# #     # return(mass_candi)
# #     if mass_candi ==[]:
# #         return(precursormz)
# #     else:
# #         max_index = intensity_candi.index(max(intensity_candi))
# #         return(mass_candi[max_index])
#
#     # return (mass_candi, max_index)
# # def get_parentpeak(mass, intensity, parention, mass_error = 0.01, ifppm = False):
# #     parention = float(parention)
# #     tolerance = so.set_tolerance(mass_error, ifppm, parention)
# #     lower_bound = parention-tolerance
# #     upper_bound = parention+tolerance
#
# #     lower_bound_i = bisect.bisect_left(mass, lower_bound)
# #     upper_bound_i = bisect.bisect_right(mass, upper_bound, lo=lower_bound_i)
#
# #     mass_candi = mass[lower_bound_i:upper_bound_i]
# #     # print("checkpoint!")
# #     intensity_candi = intensity[lower_bound_i:upper_bound_i]
# #     # return(mass_candi)
# #     if mass_candi !=[]:
# #         max_index = intensity_candi.index(max(intensity_candi))
# #         return(mass_candi[max_index], intensity_candi[max_index])
# #     else:
# #         return(parention, 0)
# def parse_formula_dict(f):
#     f = f.composition()
#     element =[]
#     count = []
#     for i in range(len(f)):
#         element.append(f[i][0])
#         count.append(f[i][1])
#     formula_dict = dict(zip(element, count))
#     return(formula_dict)

#
# # def evaluate_nl_blacklist(nl):
# #     if len(nl)==0:
# #         return(False)
# #     if len(Formula(nl).composition())==1:
# #         if Formula(nl).composition()[0][0] =='C' or Formula(nl).composition()[0][0] =='N':
# #             return(False)
# #
# #         else:
# #             return(True)
# #
# #     else:
# #         #
# #         # if:
# #         #     pass
# #         # else:
# #         return(True)
#
#
#
# # def evaluate_frag(instance):
# #     allowed_formula = []
# #     precursormz = float(spec['PrecursorMZ'])
# #     mass, intensity = prep_msms(spec['spectrum'], ifnist=True)
# #     parention = find_parention(mass, intensity, precursormz)
# #     mass, intensity = truncate_msms(mass, intensity, parention)
# #     for i in range(0, len(mass)):
# #         allowed_formula_temp = mtf.nl_to_formula(parention - mass[i], 10, spec['Formula'])
# #         if allowed_formula_temp!=[] and allowed_formula_temp!=['']:
# #             allowed_formula.extend(allowed_formula_temp)
# #             # return(allowed_formula)
# #         # else:
# #         #     allowed_formula.extend([-1])
# #             # return(allowed_formula)
# #     return (allowed_formula)
# def find_losses_nist(spec):
#     allowed_formula = []
#     precursormz = float(spec['PrecursorMZ'])
#     mass, intensity = prep_msms(spec['spectrum'], ifnist=True)
#     parention = get_parentpeak(mass, intensity, precursormz)
#     mass, intensity = truncate_msms(mass, intensity, parention)
#     for i in range(0, len(mass)):
#         allowed_formula_temp = mtf.nl_to_formula(parention - mass[i], 10, spec['Formula'])
#         if allowed_formula_temp!=[] and allowed_formula_temp!=['']:
#             allowed_formula.extend(allowed_formula_temp)
#             # return(allowed_formula)
#         # else:
#         #     allowed_formula.extend([-1])
#             # return(allowed_formula)
#     return (allowed_formula)
#
#
#
# def denoise_dyn(msms, precursor_mz):
#     msms = so.sort_spectrum(msms)
#     mass, intensity = so.break_spectra(msms)
#     index_start = np.searchsorted(mass, precursor_mz-1.6,side = 'left')
#     mass_parent = mass[index_start:]
#     intensity_parent = intensity[index_start:]
#     mass_frag = mass[0:index_start]
#     intensity_frag = intensity[0:index_start]
#     threshod =np.quantile(intensity_frag, 0.1)
#     mass_d = []
#     intensity_d = []
#     for i in range(0, len(mass_frag)):
#         if intensity_frag[i]>threshod:
#             mass_d.append(mass_frag[i])
#             intensity_d.append(intensity_frag[i])
#     if len(mass_d)==0:
#         return np.NAN
#     mass_d.extend(mass_parent)
#     intensity_d.extend(intensity_parent)
#     return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)))
# def denoise_bp(msms, precursor_mz,threshold = 1):
#     msms = so.sort_spectrum(msms)
#     mass, intensity = so.break_spectra(msms)
#     index_start = np.searchsorted(mass, precursor_mz-1.6,side = 'left')
#     mass_parent = mass[index_start:]
#     intensity_parent = intensity[index_start:]
#     mass_frag = mass[0:index_start]
#     intensity_frag = intensity[0:index_start]
#     bp_int = np.max(intensity_frag)*threshold/100
#     mass_d = []
#     intensity_d = []
#     for i in range(0, len(mass_frag)):
#         if intensity_frag[i]>bp_int:
#             mass_d.append(mass_frag[i])
#             intensity_d.append(intensity_frag[i])
#     if len(mass_d)==0:
#         return np.NAN
#     mass_d.extend(mass_parent)
#     intensity_d.extend(intensity_parent)
#     return(so.sort_spectrum(so.pack_spectra(mass_d, intensity_d)))

    # else:

    # pep_mass, pep_intensity = get_parentpeak(mass, intensity, precursormz)
    # mass, intensity = truncate_msms(mass, intensity, pep_mass)
    # if adduct=='[M+Na]+':
    #     formula = formula+'Na'
    # if adduct=='[M+NH4+]':
    #     formula = formula+'NH3'
    # if adduct=='[M+Cl]-':
    #     formula = formula+'Cl'
    # mass_d = []
    # intensity_d = []
    # threshold = so.set_tolerance(mass_error, ifppm = ifppm, precursormz=precursormz)
    #
    # # print(threshold)
    # for i in range(0,len(mass)):
    #     # print("i am evaluating", mass[i])
    #     nl_candidates = mtf.nl_to_formula(pep_mass - mass[i], threshold, formula)
    #     if nl_candidates != [] and nl_candidates!=['']:
    #         for nl in nl_candidates:
    #             if evaluate_nl_blacklist(nl)==True:
    #                 mass_d.append(mass[i])
    #                 intensity_d.append(intensity[i])
    #                 break
    # if pep_intensity != 0:
    #     mass_d.append(pep_mass)
    #     intensity_d.append(pep_intensity)
    # return(so.sort_spectra(so.pack_spectra(mass_d, intensity_d)) )