import numpy as np
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import Chem
import chemparse
import itertools
from tqdm import tqdm
import multiprocessing as mp
import re
from rdkit import RDLogger                                                                                                                                                               
RDLogger.DisableLog('rdApp.*')                                                                                                                                                           

from molmass import Formula

import spectral_denoising.spectral_operations as so
from .chem_utils import replace_adduct_string, calculate_precursormz
from spectral_denoising.constant import proton_mass
# key functions: electronic denoising and formula denoising
_numpy_formula_format = np.int16
def spectra_denoising_batch(msms_quene, smiles_quene, adduct_quene, max_allowed_deviation = 0.005):
    if len(msms_quene) != len(smiles_quene) or len(msms_quene) != len(adduct_quene):
        print('The length of msms, smiles and adduct should be the same')
        return ()
    with mp.Pool(processes=6) as pool:
            # Use starmap to handle multiple parameters
        results = pool.starmap(spectra_denoising, tqdm([(msms, smiles, adduct, max_allowed_deviation) for msms, smiles, adduct 
                                                        in zip(msms_quene, smiles_quene,adduct_quene)], total=len(msms_quene)))

    return results
def spectra_denoising(msms, smiles, adduct, max_allowed_deviation = 0.005):
    if isinstance(msms, float) or isinstance(smiles, float) or isinstance(adduct, float):
        print('the input is invalid')
        return np.nan
    
    electronic_denoised = electronic_denoising(msms)
    if isinstance(electronic_denoised, float):
        return(np.nan)
    formula_denoised = formula_denoising(electronic_denoised, smiles, adduct, max_allowed_deviation)
    return formula_denoised
def formula_denoising(msms, smiles, adduct, max_allowed_deviation=0.005):
    master_formula = prep_formula(smiles, adduct)
    msms = so.sort_spectra(msms)
    if isinstance(master_formula, float):
        print(f'Error: invalid smiles {smiles} or invalid adduct {adduct}')
        return np.nan
    computed_pmz = calculate_precursormz(adduct, smiles)
    pmz, mass_threshold = get_pmz_statistics(msms, computed_pmz, max_allowed_deviation)
    frag_msms, pmz_msms = so.slice_spectra(msms, pmz-1.6)
    all_possible_candidate_formula,all_possible_mass = get_all_subformulas(master_formula)
    denoise_tag = get_denoise_tag(frag_msms, all_possible_candidate_formula, all_possible_mass, pmz, has_benzene(smiles), mass_threshold)
    frag_msms_denoised = frag_msms[denoise_tag]
    return so.add_spectra(frag_msms_denoised, pmz_msms)

def electronic_denoising(msms):
    mass, intensity = so.break_spectra(msms)
    order = np.argsort(intensity)
    mass = mass[order]
    intensity = intensity[order]
    mass_confirmed = np.array([])
    intensity_confirmed = np.array([])
    while len(intensity)>0:
        seed_intensity = np.max(intensity)
        idx_left = np.searchsorted(intensity, seed_intensity*0.999, side= 'left')
        mass_temp = mass[idx_left:]
        intensity_temp = intensity[idx_left:]
        if len(mass_temp)<=3:
            mass_confirmed =  np.concatenate((mass_confirmed, mass_temp))
            intensity_confirmed = np.concatenate((intensity_confirmed,intensity_temp))
        intensity = intensity[0:idx_left]
        mass = mass[0:idx_left]
    if len(mass_confirmed)==0:
        return so.pack_spectra([],[])
    return(so.sort_spectra(so.pack_spectra(mass_confirmed, intensity_confirmed)) )


from .chem_utils import determine_adduct_charge, determine_parent_coefs, parse_adduct
from .seven_golden_rules import check_ratio
def get_denoise_tag(frag_msms, all_possible_candidate_formula, all_possible_mass, pmz,has_benzene, mass_threshold):
    tag = []
    if has_benzene:
        pmz = pmz+Formula('N2O').isotope.mass

    for f in frag_msms:
        loss = pmz - f[0] # get loss
        idx_left, idx_right = all_possible_mass.searchsorted([loss - mass_threshold, loss + mass_threshold])
        tag.append(check_candidates(all_possible_candidate_formula[idx_left:idx_right]))
    return tag
def check_candidates(candidates):
    for c in candidates:
        if check_ratio(c):
            return True
    return False
def get_pmz_statistics(msms, c_pmz, max_allowed_deviation):
    msms_T = msms.T
    idx_left, idx_right = msms_T[0].searchsorted([c_pmz - 0.01, c_pmz + 0.01])

    if idx_left == idx_right:
        return c_pmz, max_allowed_deviation
    pmz_idx = np.argmax(msms_T[1][idx_left:idx_right])
    r_pmz = msms_T[0][idx_left:idx_right][pmz_idx]
    r_deviation = np.abs(c_pmz - r_pmz)
    if r_deviation*1.75>max_allowed_deviation:
        return r_pmz, r_deviation*1.75
    else:
        return r_pmz, max_allowed_deviation
def get_all_subformulas(raw_formula):
    if is_smiles(raw_formula) == True:
        raw_formula = CalcMolFormula(Chem.MolFromSmiles(raw_formula))
    master_formula = chemparse.parse_formula(raw_formula)
    formula_range = [range(int(x) + 1) for (x) in master_formula.values()]
    mass_arr = [Formula(x).isotope.mass for (x) in master_formula.keys()]
    all_possible_candidate_formula_arr = np.array(list(itertools.product(*formula_range)), _numpy_formula_format)
    all_possible_mass = np.sum(mass_arr * all_possible_candidate_formula_arr, axis=1)
    order = np.argsort(all_possible_mass)
    all_possible_mass = all_possible_mass[order]
    all_possible_candidate_formula_arr = all_possible_candidate_formula_arr[order]
    all_possible_candidate_formula = [dict_to_formula(x, list(master_formula.keys())) for x in all_possible_candidate_formula_arr]
    return all_possible_candidate_formula,all_possible_mass
def dict_to_formula(candidate, element_dict):
    string = ''
    for i in range(0, len(candidate)):
        if candidate[i]>1:
            string += element_dict[i] + str(candidate[i])
        elif candidate[i]==1:
            string += element_dict[i]
    return string
def prep_formula(smiles, adduct):
    if smiles != smiles or adduct != adduct or 'i' in adduct:
        return np.nan
    adduct = replace_adduct_string(adduct)
    mol = Chem.MolFromSmiles(smiles)
    formula = CalcMolFormula(mol)
    extra_atoms = has_benzene(mol)
    if adduct in ['[M]+', '[M]-']:
        if formula[-1]==adduct[-1]:
            formula = formula[0:-1]
            master_formula = Formula(formula)
            if extra_atoms ==True:
                master_formula = master_formula.__add__(Formula('N2O'))
            return(master_formula.formula)
        else:
            print(f'the correct master formula cannot be determined for {smiles} and {adduct}')
            return(np.nan)
    charge = determine_adduct_charge(adduct)
    if abs(charge) > 1:
        print(f'the correct master formula cannot be determined for adduct charge > 1')
        return(np.nan)
    m_coef = determine_parent_coefs(adduct)
    master_formula = Formula(formula*m_coef)
    parsed_adduct = parse_adduct(adduct)
    
    for p in parsed_adduct:
        sign, count, ion_type = p
        if ion_type == 'H':
            continue # skip proton
        if sign == '+':
            master_formula = master_formula.__add__(Formula(ion_type*count))
        elif sign == '-':
            master_formula = master_formula.__sub__(Formula(ion_type*count))
        else:
            continue
    if extra_atoms == True:
        master_formula = master_formula.__add__(Formula('N2O'))

    return(master_formula.formula)
from .identifier_utils import is_mol, is_smiles
def has_benzene(molecule):
    if is_mol(molecule) == False and is_smiles(molecule) == True:
        molecule = Chem.MolFromSmiles(molecule)
    benzene = Chem.MolFromSmiles('c1ccccc1')  # Aromatic benzene SMILES notation

    # Check if benzene is a substructure of the given molecule
    return molecule.HasSubstructMatch(benzene)
# below are benchmarked functions


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

    return(so.sort_spectra(so.pack_spectra(mass_signal, intensity_signal)))
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

    return so.sort_spectra(so.pack_spectra(mass_taken, intensity_taken))
def threshold_denoising(msms, threshold = 1):
    mass, intensity = so.break_spectra(msms)
    intensity_percent = intensity/np.max(intensity)
    to_keep = intensity_percent>(threshold/100)
    mass = mass[to_keep]
    intensity = intensity[to_keep]
    return(so.pack_spectra(mass, intensity))