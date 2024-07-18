#!/bin/env/python3
import itertools
import numpy as np
import pandas as pd
import re
import math
import sys
from molmass import Formula
from toolsets.constants import single_charged_adduct_mass, atom_valence
import toolsets.denoising_related_functions as de
# from toolsets.denoising_related_functions import evaluate_nl_blacklist
from brainpy import isotopic_variants
from toolsets.search import string_search
import toolsets.spectra_operations as so
# from toolsets.spectra_operations import normalize_spectrum, break_spectra, pack_spectra, set_tolerance, standardize_spectra, num_peaks
from sklearn.linear_model import LinearRegression
atom_mass = {
    "H": 1.007825,#0
    "C": 12.00000,#1
    "N": 14.003074,#2
    "O": 15.994915,#3
    'F': 18.99840322,#4
    'Si': 27.97693,
    "P": 30.973762,#6
    "S": 31.972071,#7
    "Cl": 34.968852721,  # Cl 37: 36.9659#8
    "Br": 78.9183,  # Br 81: 80.9163#10
    # 'I': 126.9045#11

    # 'Co': 58.93319429
}
atom_list = [
    'H', 'C', 'N', 'O', 'F','Si','P', 'S', 'Cl', 'Br'
             ]
atom_valence = {
    "H": 1,
    "C": 4,
    "N": 5,
    "O": 6,
    'F': 7,
    'Si':4,
    "P": 5,
    "S": 6,
    "Cl": 7,  # Cl 37: 36.9659
    "Br": 7,  # Br 81: 80.9163
    # 'I':7
    # 'Na':22.989769282
    # 'Co': 58.93319429
}

numpy_formula_format = np.int16

atom_dict = {a: i for i, a in enumerate(atom_list)}
len_atom_dict = len(atom_dict)
atom_mass_array = np.zeros(len_atom_dict, np.float32)
for atom in atom_dict:
    atom_mass_array[atom_dict[atom]] = atom_mass[atom]


class MolecularFormula(object):
    __slots__ = ['_data', '_hash']

    def __init__(self, data=None):
        if data is not None:
            self._data = np.array(
                data, numpy_formula_format, copy=True)
        else:
            self._data = np.zeros(len_atom_dict, numpy_formula_format)
        self._hash = None
        pass

    def __getitem__(self, item):
        return self._data[atom_dict[item]]

    def __setitem__(self, key, value):
        self._data[atom_dict[key]] = value

    def __str__(self):
        string = ''

        for atom in atom_list:
            atom_num = self[atom]
            if atom_num:
                if atom_num > 1:
                    string += atom + str(atom_num)
                else:
                    string += atom
        return string

    def from_string(self, key_string):
        f = Formula(key_string)
        # f.composition()
        all_atom_nums = []
        for i in range(len(f.composition())):
            all_atom_nums.append((f.composition()[i][0], f.composition()[i][1]))
        for atom_num in all_atom_nums:
            self[atom_num[0]] = int(atom_num[1])
        pass

    def get_data(self):
        return self._data

    def get_mass(self):
        return np.sum(atom_mass_array * self._data)
def get_atom_number(formula, atom):
    if type(formula) is not MolecularFormula:
        mol = MolecularFormula()
        mol.from_string(formula)
    else:
        mol = formula
    if atom in atom_list:
        n_atom = mol.get_data()[atom_list.index(atom)]
        return(n_atom)
    else:
        print('the element is not in list!')
        return np.NAN
#
def check_chemistry(formula):
    mol = MolecularFormula()
    mol.from_string(formula)
    if _check_element_numbers(formula) and _check_senior(formula) and _check_HC_ratio(formula) and _check_hetero_ratio(formula) and _check_heuristic(formula):
        return(True)
    else:
        return False
def _check_element_numbers(formula):
    mol = MolecularFormula()
    mol.from_string(formula)
    mol_mass = mol.get_mass()
    element_upper = [
        [39,72,20,20,9,10,16,10,4,8],
        [78,126,20,27,9,14,34,12,8,14],
        [156,180,20,40,9,14,48,12,10,15]
    ]
    if mol_mass<500:
        rule = element_upper[0]
    elif mol_mass>=500 and mol_mass<=1000:
        rule = element_upper[1]
    else:
        rule = element_upper[2]
    rule_1_atoms = ['C','H','N','O','P', 'S','F','Cl', 'Br', 'Si']
    for idx in range(0, len(rule_1_atoms)):
        if get_atom_number(mol, rule_1_atoms[idx])>rule[idx]:
            return(False)
    return(True)
def _check_senior(formula):
    state = True
    mol = MolecularFormula()
    mol.from_string(formula)
    max_valence = 0
    sum_valence = 0
    sum_atom = mol.get_data().sum()
    for atom in atom_list:
        sum_valence=sum_valence+get_atom_number(mol, atom)*atom_valence[atom]
        # print(max_valence)
        if get_atom_number(mol, atom)>0:

            if atom_valence[atom] >max_valence:
                max_valence=atom_valence[atom]
    # print(sum_atom, sum_valence)
    if sum_valence%2==0:
        state = True
    elif sum_valence%2==1 and sum_atom%2==0:
        state = True
    else:
        # print('hiii')
        return(False)
    if sum_valence<2*max_valence:
        return(False)
    if sum_valence<2*sum_atom-1:
        # print('hii')
        return(False)

    return True
def _check_HC_ratio(formula):
    n_C = get_atom_number(formula, 'C')
    n_H = get_atom_number(formula, 'H')
    ratio=(n_H/n_C)
    if ratio>3.1:
        return(False)
    else:
        return (True)
def _check_heuristic(formula):

    n_N = get_atom_number(formula, 'N')
    n_O = get_atom_number(formula, 'O')
    n_P = get_atom_number(formula, 'P')
    n_S = get_atom_number(formula, 'S')
    if n_N>1 and n_O >1 and n_P>1 and n_S >1:
        if n_N >10 or n_O>20 or n_P >4 or n_S>3:
            return(False)
    if n_N >3 and n_O>3 and n_P >3:
        if n_N >11 or n_O >22 or n_P >6:
            return(False)
    if n_O >1 and n_P >1 and n_S >1:
        if n_O >14 or n_P >3 or n_S >3:
            return(False)
    if n_P >1 and n_S >1 and n_N >1:
        if n_P >3 or n_S >3 or n_N >4:
            return(False)
    if n_N >6 and n_O >6 and n_S>6:
        if n_N >19 or n_O >14 or n_S>8:
            return False
    return True
def _check_hetero_ratio(formula):
    n_C = get_atom_number(formula, 'C')
    n_H = get_atom_number(formula, 'H')
    n_F = get_atom_number(formula, 'F')
    n_Cl = get_atom_number(formula, 'Cl')
    n_Br = get_atom_number(formula, 'Br')
    n_N = get_atom_number(formula, 'N')
    n_O = get_atom_number(formula, 'O')
    n_P = get_atom_number(formula, 'P')
    n_S = get_atom_number(formula, 'S')
    n_Si = get_atom_number(formula, 'Si')
    if n_F/n_C >6:
        return (False)
    if n_Cl/n_C >0.8:
        return False
    if n_Br>n_C >0.8:
        return False
    if n_N/n_C>1.3:
        return False
    if n_O/n_C>1.2:
        return False
    if n_P/n_C >0.3:
        return False
    if n_S/n_C >0.8:
        return False
    if n_Si/n_C >0.5:
        return False
    return True


def precursor_to_formula(precursor_mz, adduct, mass_error = 0.005):
    base_mass = precursor_mz-single_charged_adduct_mass[adduct]
    all_candidates = mass_to_formula(base_mass, mass_error)
    return(all_candidates)
def standardize_formula(formula):
    t = Formula(formula)

    return (t.formula)
def mass_to_formula(mass, mass_error=0.005):
    mol_mass = mass
    lo_mass = mol_mass - mass_error
    hi_mass = mol_mass + mass_error
    result = []
    candidate_formula_array = np.zeros(len_atom_dict, numpy_formula_format)
    _calculate_formula(lo_mass, hi_mass, candidate_formula_array,
                       len(candidate_formula_array) - 1, result)
    formulas = []
    for r in result:
        formulas.append(standardize_formula(r.__str__()))

    return formulas
def _calculate_formula(mass_start, mass_end, candidate_formula_array, cur_i, result):
    # print(cur_i)
    atom_mass_cur = atom_mass_array[cur_i]
    # print(atom_mass_cur)
    atom_num = math.floor(mass_end / atom_mass_cur)
    if cur_i == 0:
        # print('i am in if loop')
        # This is H
        h_num_low = mass_start / atom_mass_cur
        if atom_num >= h_num_low:
            # print("i am in ifif")
            candidate_formula_array[0] = atom_num
            result.append(MolecularFormula(candidate_formula_array))
        # else:
        #     candidate_formula_array[0] = atom_num-1
        #     result.append(MolecularFormula(candidate_formula_array))
    else:
        # print('i am in else loop')
        if atom_num>0:
            for i in range(atom_num+1):
                f = np.copy(candidate_formula_array)
                f[cur_i] = i
                # print(mass_start - i * atom_mass_cur)
                # print(atom_list[cur_i])
                _calculate_formula(mass_start - i * atom_mass_cur, mass_end - i * atom_mass_cur,
                                   f, cur_i - 1, result)

        else:
            f = np.copy(candidate_formula_array)
            f[cur_i] = 0
            _calculate_formula(mass_start, mass_end,
                               f, cur_i - 1, result)
# # M, f -> measured, m, p-> therotical
def ms1_intensity_score(ms1_silico, ms1, iso_state = 2):
    # ms1_silico =standardize_spectra(ms1_silico)
    # ms1 = standardize_spectra(ms1)
    from toolsets.spectra_operations import break_spectra
    mass_silico, intensity_silico = break_spectra(ms1_silico)
    mass, intensity = break_spectra(ms1)
    score= 1

    for cur_i in range(0, iso_state+1):
        if abs(mass_silico[cur_i]-mass[cur_i])>0.5:
            return 0
        score = score*(1-abs(intensity[cur_i]-intensity_silico[cur_i]))
    return score
def ms1_intensity_score(ms1_silico, ms1, iso_state = 2):
    # ms1_silico =standardize_spectra(ms1_silico)
    # ms1 = standardize_spectra(ms1)
    from toolsets.spectra_operations import break_spectra
    mass_silico, intensity_silico = break_spectra(ms1_silico)
    mass, intensity = break_spectra(ms1)
    score= 1

    for cur_i in range(0, iso_state+1):
        if abs(mass_silico[cur_i]-mass[cur_i])>0.5:
            return 0
        score = score*(1-abs(intensity[cur_i]-intensity_silico[cur_i]))
    return score
def bin_ms1(ms1_temp, precursor_mz):
    ms1_temp = so.sort_spectrum(ms1_temp)
    mass,intensity = so.break_spectra(ms1_temp)
    offset = [abs(precursor_mz-x) for x in mass]
    if np.min(offset)>0.5:
        return(np.NAN)
    m0_idx = np.argmin(offset)
    m0 = find_isotope(ms1_temp, precursor_mz, 0)
    m1 = find_isotope(ms1_temp, precursor_mz, 1)
    m2 = find_isotope(ms1_temp, precursor_mz, 2)
    intensity_binned = []
    mass_binned = []
    for i in [m0, m1, m2]:
        idx_lower = np.searchsorted(mass, i[0]-0.02, side='left')
        idx_upper = np.searchsorted(mass, i[0]+0.02, side = 'right')
        intensity_binned.append(np.sum(intensity[idx_lower:idx_upper]))
        mass_binned.append(i[0])

    ms1_binned = so.pack_spectra(mass_binned,intensity_binned)
    ms1_binned = so.standardize_spectra(ms1_binned)
    return(ms1_binned)



def find_isotope(ms1_temp, precursor_mz, state):
    mass, intensity = so.break_spectra(ms1_temp)
    mn = precursor_mz+state
    offset = [abs(mn -x) for x in mass]
    if np.min(offset)>0.5:
        return(mn, 0)
    else:
        return(mass[np.argmin(offset)], intensity[np.argmin(offset)])
#
# def get_spec_p_intensity(ms1_silico, ms1, off_set = 0.002):
#     mj, pj = break_spectra(ms1_silico)
#     Mj,  fj = break_spectra(ms1)
#     fj = [item +off_set for item in fj]
#     fj = [item/sum(fj) for item in fj]
#     p_intensity = 1
#     for cur_i in range(0, len(fj)):
#         # print(get_p_intensity(fj[cur_i], pj[cur_i]))
#         p_intensity = p_intensity*get_p_intensity(fj[cur_i], pj[cur_i])
#         # print(p_intensity)
#     return(p_intensity)
# def get_spec_p_mass(ms1_silico, ms1, include_base = False):
#     mj, pj = break_spectra(ms1_silico)
#     Mj,  fj = break_spectra(ms1)
#     M0 = Mj[0]
#     m0 = mj[0]
#     p_mass = 1
#     if include_base ==True:
#         start = 0
#     else:
#         start = 1
#     for cur_i in range(start, len(Mj)):
#         if cur_i ==0:
#             p_mass = p_mass*get_p_mass(Mj[cur_i], mj[cur_i], fj[cur_i])
#         else:
#             p_mass = p_mass*get_p_mass(Mj[cur_i], mj[cur_i], fj[cur_i], M0 = M0, m0 = m0)
#     return(p_mass)
# def get_p_intensity(fj, pj):
#     # fj is the intensity of measured peak
#     sigma = get_sigma_intensity(fj)
#     p = math.erfc(abs(math.log( fj/pj, 10))/(math.sqrt(2)*sigma))
#     return(p)
#
# def get_p_mass(Mj, mj, fj, M0 = None, m0 = None):
#     sigma = get_sigma_mass(fj)
#     # print(sigma)
#     # print('alpha1, alpha 0',alpha1, alpha0)
#     if M0 !=None and m0!=None:
#         # calculate by difference to eliminate sysmtematic error
#         Mj = Mj-M0
#         mj = mj-m0
#     #     this is when comaring monoisotopic peaks
#     p = math.erfc(abs(Mj-mj)/(math.sqrt(2)*sigma))
#     return(p)
# def get_sigma_mass(intensity, mass_accuray = 0.005):
#     alpha = get_alpha(intensity)
#     return(1/3*alpha*mass_accuray)
# def get_sigma_intensity(intensity):
#     beta = get_beta(intensity)
#     return(1/3*math.log((1+beta),10 ))
#
# def get_beta(intensity):
#     if intensity >1:
#         intensity = intensity/100
#     if intensity >=0.8:
#         return (0.08)
#     elif intensity >=0.5 and intensity < 0.8:
#         x = np.array([0.5, 0.8]).reshape((-1, 1))
#         y = np.array([0.12, 0.08])
#
#     elif intensity >=0.2 and intensity <0.5:
#         x = np.array([0.2, 0.5]).reshape((-1, 1))
#         y = np.array([0.25, 0.12])
#     elif intensity >=0.1 and intensity <0.2:
#         x = np.array([0.1,0.2 ]).reshape((-1, 1))
#         y = np.array([0.5, 0.25])
#     elif intensity >=0.05 and intensity <0.1:
#         x = np.array([0.05, 0.1]).reshape((-1, 1))
#         y = np.array([1.2, 0.5])
#     elif intensity >=0.01 and intensity <0.05:
#         x = np.array([0.01, 0.05]).reshape((-1, 1))
#         y = np.array([2, 1.2])
#     elif intensity <0.01:
#         return(2)
#     model = LinearRegression().fit(x, y)
#     return (model.predict(np.array([intensity]).reshape(1,-1))[0])
# def get_alpha(intensity):
#     if intensity >1:
#         intensity = intensity/100
#     if intensity >=0.8:
#         return (1)
#     elif intensity >=0.5 and intensity < 0.8:
#         x = np.array([0.8, 0.5]).reshape((-1, 1))
#         y = np.array([1, 1.1])
#     elif intensity >=0.2 and intensity <0.5:
#         x = np.array([0.5, 0.2]).reshape((-1, 1))
#         y = np.array([1.1, 1.2])
#     elif intensity >=0.1 and intensity <0.2:
#         x = np.array([0.2, 0.1]).reshape((-1, 1))
#         y = np.array([1.2, 2])
#     elif intensity >=0.05 and intensity <0.1:
#         x = np.array([0.05, 0.1]).reshape((-1, 1))
#         y = np.array([2.5, 2])
#     elif intensity >=0.01 and intensity <0.05:
#         x = np.array([0.01, 0.05]).reshape((-1, 1))
#         y = np.array([2.5, 10])
#     elif intensity <0.01:
#         return(10)
#     model = LinearRegression().fit(x, y)
#     return (model.predict(np.array([intensity]).reshape(1,-1))[0])
#
# def get_candidate(mass, mass_error, ifppm):
#     result = neutral_to_formula(mass,
#         set_tolerance(mass_error=mass_error,  ifppm=ifppm, precursormz=mass))
#     compliance = []
#     formulas = []
#     for i in (range(len(result))):
#         compliance.append(check_chemistry(result[i].__str__()))
#         # f_temp = Formula
#         formulas.append(Formula(result[i].__str__()).formula)
#         # RDBE.append(get_RDBE(Formula(result[i].__str__())))
#     tentative_result = pd.DataFrame(list(zip(formulas, compliance )), columns=['Formulas', 'chemistry'])
#     tentative_result = string_search(tentative_result, 'chemistry', True)
#     tentative_result.reset_index(inplace = True, drop = True)
#     return(tentative_result)
# atom_dict = {a: i for i, a in enumerate(atom_list)}
# len_atom_dict = len(atom_dict)
# atom_mass_array = np.zeros(len_atom_dict, np.float32)
# for atom in atom_dict:
#     atom_mass_array[atom_dict[atom]] = atom_mass[atom]
#
#
# def neutral_to_formula(mass, mass_error):
#     mass = float(mass)
#     mass_error = float(mass_error)
#     return precursor_mass_to_formula(mass+1.007276, mass_error, '[M+H]+')

def parse_formula(f):
    f = f.composition()
    element =[]
    count = []
    for i in range(len(f)):
        element.append(f[i][0])
        count.append(f[i][1])
    df = pd.DataFrame(list(zip(element, count)),
               columns =['element', 'count'])
    return(df)
# def get_nvalence(formula):
#     # this formula is already a molmass.Formula item!!!!
#     df = parse_formula(formula)
#     total_valence = 0
#     # return(df)
#     for index, row in df.iterrows():
#         total_valence=total_valence+atom_valence[row['element']]*row['count']
#     return(total_valence)
#
# def get_num_atom(df, atoms):
#
#     # atoms has to be in list
#     df_atoms = df[df['element'].isin(atoms)]
#     return(df_atoms['count'].sum())
#
# def get_RDBE(formula):
#     # formula = Formula(formula)
#     df = parse_formula(formula)
#     RDBE = get_num_atom(df, ['C'])-(get_num_atom(df, ['H'])+get_num_atom(df, ['F','Cl','Br','I']))/2+get_num_atom(df, ['N'])/2+1
#     return RDBE
#
# def check_senior3(formula):
#     # print("i am in new")
#
#     # f = Formula(formula)
#     # RDBE = get_RDBE(f)
#     natoms = formula.atoms
#     nvalence = get_nvalence(formula)
#     if nvalence >= 2*natoms-1:
#         return (True)
#     else:
#         return (False)
#
# def check_chemistry(formula):
#     formula = Formula(formula)
#     nvalence = get_nvalence(formula)
#     if check_senior3(formula) == True and get_RDBE(formula) >-0.5 and de.evaluate_nl_blacklist(formula.formula) and nvalence%2==0:
#         return (True)
#     else:
#         return (False)
#
from brainpy import isotopic_variants
from toolsets.constants import single_charged_adduct_mass
def get_isotope(formula, digit = 3, adduct = None):
    from toolsets.spectra_operations import pack_spectra, standardize_spectra
    formula_temp = Formula(formula)
    formula_temp = parse_formula(formula_temp)
    mol_temp = dict(zip(formula_temp['element'], formula_temp['count']))
    theoretical_isotopic_cluster = isotopic_variants(mol_temp, npeaks=digit, charge=0)
    mass = []
    intensity = []
    for i in range(digit):
        if adduct != None:
            mass.append(theoretical_isotopic_cluster[i].mz+single_charged_adduct_mass[adduct])
        else:
            mass.append(theoretical_isotopic_cluster[i].mz)
        intensity.append(theoretical_isotopic_cluster[i].intensity)
    ms1 = pack_spectra(mass, intensity)
    ms1 = standardize_spectra(ms1)
    return (ms1)
#
# def _calculate_formula(mass_start, mass_end, candidate_formula_array, cur_i, result):
#     # print(cur_i)
#     atom_mass_cur = atom_mass_array[cur_i]
#     # print(atom_mass_cur)
#     atom_num = math.floor(mass_end / atom_mass_cur)
#     if cur_i == 0:
#         # print('i am in if loop')
#         # This is H
#         h_num_low = mass_start / atom_mass_cur
#         if atom_num >= h_num_low:
#             # print("i am in ifif")
#             candidate_formula_array[0] = atom_num
#             result.append(MolecularFormula(candidate_formula_array))
#         # else:
#         #     candidate_formula_array[0] = atom_num-1
#         #     result.append(MolecularFormula(candidate_formula_array))
#     else:
#         # print('i am in else loop')
#         if atom_num>0:
#             for i in range(atom_num+1):
#                 f = np.copy(candidate_formula_array)
#                 f[cur_i] = i
#                 # print(mass_start - i * atom_mass_cur)
#                 # print(atom_list[cur_i])
#                 _calculate_formula(mass_start - i * atom_mass_cur, mass_end - i * atom_mass_cur,
#                                    f, cur_i - 1, result)
#
#         else:
#             f = np.copy(candidate_formula_array)
#             f[cur_i] = 0
#             _calculate_formula(mass_start, mass_end,
#                                    f, cur_i - 1, result)
#
# def precursor_mass_to_formula(mass, mass_error, addition):
#     mol_mass = mass - single_charged_adduct_mass[addition]
#     # print(mol_mass)
#     lo_mass = mol_mass - mass_error
#     hi_mass = mol_mass + mass_error
#     # print(mass_error)
#
#     result = []
#     candidate_formula_array = np.zeros(len_atom_dict, numpy_formula_format)
#     _calculate_formula(lo_mass, hi_mass, candidate_formula_array,
#                        len(candidate_formula_array) - 1, result)
#     return result
#
#
#
#
#
#
# class MolecularFormula(object):
#     __slots__ = ['_data', '_hash']
#
#     def __init__(self, data=None):
#         if data is not None:
#             self._data = np.array(
#                 data, numpy_formula_format, copy=True)
#         else:
#             self._data = np.zeros(len_atom_dict, numpy_formula_format)
#         self._hash = None
#         pass
#
#     def __getitem__(self, item):
#         return self._data[atom_dict[item]]
#
#     def __setitem__(self, key, value):
#         self._data[atom_dict[key]] = value
#
#     def __str__(self):
#         string = ''
#
#         for atom in atom_list:
#             atom_num = self[atom]
#             if atom_num:
#                 if atom_num > 1:
#                     string += atom + str(atom_num)
#                 else:
#                     string += atom
#         return string
#
#     def from_string(self, key_string):
#         f = Formula(key_string)
#         # f.composition()
#         all_atom_nums = []
#         for i in range(len(f.composition())):
#             all_atom_nums.append((f.composition()[i][0], f.composition()[i][1]))
#         for atom_num in all_atom_nums:
#             self[atom_num[0]] = int(atom_num[1])
#         pass
#
#     def get_data(self):
#         return self._data
#
#     def get_mass(self):
#         return np.sum(atom_mass_array * self._data)

# def product_mass_to_formula(mass, precursor_formula, mass_error = 0.01, adduct=None ):
#     if adduct != None:
#         mass -= single_charged_adduct_mass[adduct]

#     lo_mass = mass - mass_error
#     hi_mass = mass + mass_error

#     # Generate candidate range
#     precursor_data = MolecularFormula()
#     precursor_data.from_string(precursor_formula)
#     formula_range = [range(x + 1) for x in precursor_data.get_data()]
#     all_possible_candidate_formula = np.array(
#         list(itertools.product(*formula_range)), numpy_formula_format)
#     all_possible_mass = np.sum(
#         atom_mass_array * all_possible_candidate_formula, axis=1)

#     candidate_data = all_possible_candidate_formula[(lo_mass <= all_possible_mass) & (all_possible_mass <= hi_mass)]

#     result = []
#     for data in candidate_data:
#         formula = MolecularFormula(data)
#         result.append(formula)
#     return result