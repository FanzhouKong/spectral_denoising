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
from .identifier_utils import is_mol, is_smiles, is_formula
from molmass import Formula
from . import spectral_operations as so
from .chem_utils import replace_adduct_string, calculate_precursormz
from .constant import proton_mass
# key functions: electronic denoising and formula denoising
_numpy_formula_format = np.int16
import warnings
warnings.filterwarnings("ignore")
def spectral_denoising_batch(msms_query, smiles_query, adduct_query, mass_tolerance = 0.005):
    """
    Perform batch spectral denoising on multiple sets of MS/MS spectra, SMILES strings, and adducts. Uses multiprocessing to parallelize the denoising process.

    Parameters:
        msms_query (list): A list of MS/MS spectra data.

        smiles_query (list): A list of SMILES strings corresponding to the MS/MS spectra.

        adduct_query (list): A list of adducts corresponding to the MS/MS spectra.

        mass_tolerance (float, optional): The allowed deviation for the denoising process. Default is 0.005.
    Returns:
        list: A list of denoised MS/MS from the spectral denoising process.
    Notes:
        - The lengths of msms_query, smiles_query, and adduct_query must be the same. If not, the function will print an error message and return an empty tuple.
        - The function uses multiprocessing to parallelize the denoising process, utilizing 6 processes.
    """
    print('i am in new')
    if len(msms_query) != len(smiles_query) or len(msms_query) != len(adduct_query):
        print('The length of msms, smiles and adduct should be the same')
        return ()
    with mp.Pool(processes=6) as pool:
            # Use starmap to handle multiple parameters
        results = pool.starmap(spectral_denoising, tqdm([(msms, smiles, adduct, mass_tolerance) for msms, smiles, adduct 
                                                        in zip(msms_query, smiles_query,adduct_query)], total=len(adduct_query)))

    return results
def spectral_denoising_with_master_formulas(msms, master_formula, benzene_tag, query_pmz,mass_tolerance = 0.005):
    if isinstance(msms, float):
        return np.nan
    if isinstance(master_formula, float):
        return msms
    # msms = so.sort_spectrum(msms)
    if master_formula != master_formula:
        # print(f'Error: invalid smiles {smiles} or invalid adduct {adduct}')
        return msms
    pmz, real_mass_threshold = get_pmz_statistics(msms, query_pmz, mass_tolerance)
    mass_threshold=mass_tolerance
    # if real_mass_threshold>mass_tolerance:
    #     mass_threshold=real_mass_threshold*1.5
    # else:
    #     mass_threshold = mass_tolerance
    frag_msms, pmz_msms = so.slice_spectrum(msms, pmz-1.6)
    all_possible_candidate_formula,all_possible_mass = get_all_subformulas(master_formula)
    denoise_tag = get_denoise_tag(frag_msms, all_possible_candidate_formula, all_possible_mass, pmz, benzene_tag, mass_threshold)
    frag_msms_denoised = frag_msms[denoise_tag]
    return so.add_spectra(frag_msms_denoised, pmz_msms)
def spectral_denoising(msms, smiles, adduct, mass_tolerance = 0.005):
    """
    Perform spectral denoising on the given mass spectrometry data. The function first performs electronic denoising, followed by formula denoising.

    Parameters:
        msms (numpy.array): The mass spectrometry data to be denoised.
        smiles (str): The SMILES representation of the molecule.
        adduct (str): The adduct type.
        mass_tolerance (float, optional): The mass tolerance for the denoising process. Default is 0.005.
    Returns:
        numpy.array: The denoised mass spectrometry data.Returns NaN if the input is invalid or if the denoising process fails.
    Notes:
        - The function first checks if any of the inputs are of type np.nan, which is considered invalid.
        - It then performs electronic denoising on the msms data.
        - If electronic denoising resulted in empty spectrum (all ions removed), it will return np.nan.
        - If successful, it proceeds to formula denoising using the electronic denoised data, smiles, adduct, and mass_tolerance.
    """

    if isinstance(msms, float) or isinstance(smiles, float) or isinstance(adduct, float):
        # print('the input is invalid')
        return np.nan
    
    electronic_denoised = electronic_denoising(msms)
    if isinstance(electronic_denoised, float):
        return(np.nan)
    formula_denoised = formula_denoising(electronic_denoised, smiles, adduct, mass_tolerance)
    return formula_denoised
def formula_denoising(msms, smiles, adduct, mass_tolerance=0.005):
    """
    Perform formula denoising on the given mass spectrometry data. The function first re-generate formula based on chemical rules, get the statistic of the precursor m/z, and then perform formula denoising.
    The precursor region is not affected by the denoising process, only the frgamnet region is denoised.

    Parameters:
        msms (numpy.array): The mass spectrometry data to be denoised. 
        smiles (str): The SMILES string representing the molecular structure. This will also recognize the molecular formula as input, but risking leading to false positives due to not incorporting possibilities of forming extra N2/H2O adducts.
        adduct (str): The adduct type used in the mass spectrometry.
        mass_tolerance (float, optional): The mass tolerance for precursor m/z calculation. Default is 0.005.
    Returns:
        numpy.ndarray: The denoised mass spectrometry data, or np.nan If the SMILES string or adduct is invalid, or all ions removed.
    """
    
    master_formula = prep_formula(smiles, adduct)#check
    msms = so.sort_spectrum(msms)
    if master_formula != master_formula:
        # print(f'Error: invalid smiles {smiles} or invalid adduct {adduct}')
        return msms
    computed_pmz = calculate_precursormz(adduct, smiles)#check
    
    if computed_pmz != computed_pmz:
        return msms
    pmz, real_mass_threshold = get_pmz_statistics(msms, computed_pmz, mass_tolerance)
    if real_mass_threshold>mass_tolerance:
        mass_threshold=real_mass_threshold*1.5
    else:
        mass_threshold = mass_tolerance
    frag_msms, pmz_msms = so.slice_spectrum(msms, pmz-1.6)
    all_possible_candidate_formula,all_possible_mass = get_all_subformulas(master_formula)
    if is_smiles(smiles) == True:
        benzene_tag = has_benzene(smiles)
    else:
        benzene_tag = True
    denoise_tag = get_denoise_tag(frag_msms, all_possible_candidate_formula, all_possible_mass, pmz,benzene_tag, mass_threshold)
    frag_msms_denoised = frag_msms[denoise_tag]
    return so.add_spectra(frag_msms_denoised, pmz_msms)
def prep_formula(smiles, adduct):
    """
    Prepares the molecular formula based on the given SMILES string and adduct.

    Args:
        smiles (str): The SMILES representation of the molecule.
        adduct (str): The adduct string representing the ionization state.
    Returns:
        str: The calculated molecular formula, or NaN if the formula cannot be determined.
    """

    if smiles != smiles or adduct != adduct or 'i' in adduct:
        return np.nan
    adduct = replace_adduct_string(adduct)
    if is_smiles(smiles) == True:
        mol = Chem.MolFromSmiles(smiles)
        formula = CalcMolFormula(mol)
        extra_atoms = has_benzene(mol)
        if adduct in ['[M]+', '[M]-'] and Chem.GetFormalCharge(mol)!= 0:
            # print('i am in loop')
            if ((Chem.GetFormalCharge(mol))>0 and adduct[-1]=='+') or ((Chem.GetFormalCharge(mol))<0 and adduct[-1]=='-'):
                # print('i am in loop')
                formula = formula[0:-1]
                master_formula = Formula(formula)
                if extra_atoms ==True:
                    master_formula = master_formula.__add__(Formula('N2O'))
                    master_formula = master_formula.formula
                    return(master_formula)
                else:
                    # print(f'the correct master formula cannot be determined for {smiles} and {adduct}')
                    return(np.nan)
            elif adduct in ['[M]+', '[M]-'] and Chem.GetFormalCharge(mol)== 0:
                return np.nan
            elif Chem.GetFormalCharge(mol)!= 0 and  adduct not in ['[M]+', '[M]-']:
                return np.nan
            else:
                return np.nan
    else:#if smiles is not a smiles, then it is a formula
        formula = smiles
        extra_atoms = True #default to true, but could lead to false negatives
        # print(formula)
        if is_formula(formula) == False:
            return np.nan
        elif abs(Formula(formula).charge) >1:
            return np.nan
        
        elif (Formula(formula).charge == 1 and adduct == '[M]+') or (Formula(formula).charge == -1 and adduct == '[M]-'):
            formula = formula[0:-1]
        elif adduct in ['[M]+', '[M]-'] and Formula(formula).charge == 0:
            return np.nan
        elif Formula(formula).charge != 0 and adduct not in ['[M]+', '[M]-']:
            return np.nan
        # else:
        #     formula = formula
    charge = determine_adduct_charge(adduct)
    if abs(charge) > 1:
        # print(f'the correct master formula cannot be determined for adduct charge > 1')
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
def electronic_denoising(msms):
    """
    Perform electronic denoising on a given mass spectrometry (MS/MS) spectrum.
    This function processes the input MS/MS spectrum by sorting the peaks based on their intensity,
    and then iteratively selects and confirms peaks based on a specified intensity threshold.
    The confirmed peaks are then packed and sorted before being returned.

    Parameters:
        msms (np.ndarray): The first item is always m/z and the second item is intensity.

    Returns:
        np.ndarray: The cleaned spectrum with electronic noises removed. If no ion presents, will return np.nan.
    """
    if isinstance(msms, float):
        return np.nan
    mass, intensity = msms.T[0], msms.T[1]
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
        return np.nan
    return(so.sort_spectrum(so.pack_spectrum(mass_confirmed, intensity_confirmed)) )


from .chem_utils import determine_adduct_charge, determine_parent_coefs, parse_adduct
from .seven_golden_rules import *
def get_denoise_tag(frag_msms, all_possible_candidate_formula, all_possible_mass, pmz,has_benzene, mass_threshold):
    
    """
    Determine which ions in the fragment regions are chemically feasible ion.
    This function calculates the mass loss for each fragment in the MS/MS data and 
    searches for candidate formulas within a specified mass threshold. If the 
    `has_benzene` flag is set, the precursor mass (`pmz`) is adjusted by adding the 
    mass of the N2O isotope to count for rare cases of forming N2/H2O adducts in the collision chamber.
    The ions will be given a True only if it can be associated with at least 1 chemically feasible subformula of the molecular formula.

    Args:
        frag_msms (numpy.array): Array of fragment MS/MS data, where each tuple contains the mass and intensity of a fragment.
        all_possible_candidate_formula (list): List of all possible candidate formulas.
        all_possible_mass (numpy.ndarray): Sorted array of all possible masses.
        pmz (float): Precursor mass.
        has_benzene (bool): Flag indicating if benzene is present.
        mass_threshold (float): Mass threshold for searching candidate formulas.
    Returns:
        list: List of denoise tags for each fragment.
    """


    tag = []
    if has_benzene:
        pmz = pmz+Formula('N2O').isotope.mass

    for f in frag_msms:
        loss = pmz - f[0] # get loss
        idx_left, idx_right = all_possible_mass.searchsorted([loss - mass_threshold, loss + mass_threshold+1E-9])
        tag.append(check_candidates(all_possible_candidate_formula[idx_left:idx_right]))
    return tag
def check_candidates(candidates):
    """
    Checks a list of candidates to see if any of them meet a certain ratio condition.

    Args:
        candidates (list): A list of candidate formulas to be checked.
    Returns:
        bool: True if at least one candidate meets the ratio condition, False otherwise.
    """
    for c in candidates:
        if check_ratio(c) and check_huristic(c):
        # if check_ratio(c) and check_senior(c):
            return True
    return False
def get_pmz_statistics(msms, c_pmz, mass_tolerance):
    """
    Use the real precursor m/z to estimate the mass deviation in a given spectrum.

    Parameters:
        msms (numpy.ndarray): A 2D array where the first row contains m/z values and the second row contains intensity values.
        c_pmz (float): The computed m/z value around which to search for the most intense peak.
        mass_tolerance (float): The mass tolerance within which to search for the most intense peak.

    Returns:
        tuple: A tuple containing:
        - r_pmz (float): The actual precursor m/z. If not found (precursor is fully fragmented), the computed m/z is returned.
        - float: The deviation between computed and actual precursor m/z, scaled by 1.75 if it exceeds the initial mass tolerance.
    """

    msms_T = msms.T
    idx_left, idx_right = msms_T[0].searchsorted([c_pmz - 0.01, c_pmz + 0.01])

    if idx_left == idx_right:
        return c_pmz, mass_tolerance
    pmz_idx = np.argmax(msms_T[1][idx_left:idx_right])
    r_pmz = msms_T[0][idx_left:idx_right][pmz_idx]
    r_deviation = np.abs(c_pmz - r_pmz)
    return r_pmz, r_deviation*1.2
def get_all_subformulas(raw_formula):
    """
    Generate all possible subformulas and their corresponding masses from a given chemical formula.

    Args:
        raw_formula (str): The input chemical formula, which can be in SMILES format or a standard chemical formula.
    Returns:
        tuple: A tuple containing:
            - all_possible_candidate_formula (list of str): A list of all possible subformulas derived from the input formula.
            - all_possible_mass (numpy.ndarray): An array of masses corresponding to each subformula.
    Notes:
        - If the input formula is in SMILES format, it will be converted to a standard chemical formula.
        - The function uses the `chemparse` library to parse the chemical formula and `itertools.product` to generate all possible combinations of subformulas.
        - The resulting subformulas and their masses are sorted in ascending order of mass for enhancing search process.
    """

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
    """
    Helper function, to get the chemical formula from a candidate list and element dictionary.

    Args:
        candidate (list of int): A list where each index corresponds to an element in `element_dict` and the value at each index represents the count of that element.
        element_dict (list of str): A list of element symbols where the index corresponds to the element's position in the `candidate` list.
    Returns:
        str: A string representing the chemical formula, where each element symbol is followed by its count if greater than 1.
    """

    string = ''
    for i in range(0, len(candidate)):
        if candidate[i]>1:
            string += element_dict[i] + str(candidate[i])
        elif candidate[i]==1:
            string += element_dict[i]
    return string


def has_benzene(molecule):
    """
    Check if the given molecule contains a benzene ring.

    Args:
        molecule (Union[Chem.Mol, str]): The molecule to check. It can be a RDKit molecule object
                                         or a SMILES string.
    Returns:
        bool: True if the molecule contains a benzene ring, False otherwise.
    """

    if is_mol(molecule) == False and is_smiles(molecule) == True:
        molecule = Chem.MolFromSmiles(molecule)
    elif is_formula(molecule) == True:
        return True
    benzene = Chem.MolFromSmiles('c1ccccc1')  # Aromatic benzene SMILES notation

    # Check if benzene is a substructure of the given molecule
    return molecule.HasSubstructMatch(benzene)
# below are benchmarked functions


def dnl_denoising(msms):
    """
    Perform Dynamic noise level estimation denoising on given msms spectra.
    Details about the algorithm can be found in the paper: A Dynamic Noise Level Algorithm for Spectral Screening of Peptide MS/MS Spectra.

    Parameters:
        msms (numpy.ndarray): A 2D numpy array with shape (2, n) where n is the number of data points. For each instance, first item is pmz and second item is intensity.
    Returns:
        numpy.ndarray: A 2D numpy array containing the denoised mass spectrometry data, sorted and packed. If the input data has only two points and does not meet the criteria, returns NaN.
    Notes:
        - The function assumes that the input data is a numpy array with two columns.
        - The function uses a linear regression model to predict the signal region.

    """

    mass, intensity = msms.T[0], msms.T[1]
    order= np.argsort(intensity)
    mass = mass[order]
    intensity = intensity[order]
    from sklearn.linear_model import LinearRegression
    if intensity[1]/2>=intensity[0]*1.5:
        signal_idx = 1

    else:
        k = 2
        if len(mass)==2:
            return(np.nan)
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

    return(so.sort_spectrum(so.pack_spectrum(mass_signal, intensity_signal)))
def ms_reduce(msms, reduce_factor = 90):
    """
    Reimplementation of MS-Reduce algorithm. 
    Details about this algorithm can be found at: MS-REDUCE: an ultrafast technique for reduction of big mass spectrometry data for high-throughput processing
    
    Parameters:
        msms (numpy.ndarray): A 2D numpy array with shape (2, n) where n is the number of data points. For each instance, first item is pmz and second item is intensity.
        reduce_factor (int, optional): The percentage by which to reduce the number of peaks. Default is 90.
    Returns:
        numpy.ndarray: The reduced MS/MS spectrum as a 2D numpy array, sorted and packed.
    """

    mass, intensity = msms.T[0], msms.T[1]
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

    return so.sort_spectrum(so.pack_spectrum(mass_taken, intensity_taken))
def threshold_denoising(msms, threshold = 1):
    """
    The most widely used and simple denoising algorithm, which discard all peaks below a predefined threshold.
    This function filters out peaks in the mass spectrometry spectrum whose 
    intensity is below a specified threshold percentage of the maximum intensity.

    Parameters:
        msms (numpy.ndarray): A 2D numpy array with shape (2, n) where n is the number of data points. For each instance, first item is pmz and second item is intensity.
        threshold (float, optional): The threshold percentage (0-100) of the maximum intensity below which peaks will be removed. Default is 1.
    Returns:
        numpy.ndarray: denoised spectrum as a 2D numpy array, sorted and packed.
    """

    mass, intensity = msms.T[0], msms.T[1]
    intensity_percent = intensity/np.max(intensity)
    to_keep = intensity_percent>(threshold/100)
    mass = mass[to_keep]
    intensity = intensity[to_keep]
    return(so.pack_spectrum(mass, intensity))