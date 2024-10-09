import chemparse
import numpy as np
from molmass import Formula
from .constant import valence_dict
def check_senior(formula):
    """
    Check if a given chemical formula satisfies SENIOR rules. not used in spectral_denoising and denoising search since SENIOR rules presumes an intact molecule, while losses can
    be formed in the fragmentation process.
    This function evaluates a chemical formula to determine if it meets the following criteria:
    1. The sum of the valences of all elements must be at least twice the maximum valence of any single element.
    2. The sum of the valences of all elements must be at least twice the number of elements minus one.

    Args:
        formula (str): A string representing the chemical formula to be evaluated.
    Returns:
        bool: True if the formula satisfies the valence rules, False otherwise.
    """

    max_valence = 0
    sum_valence = 0
    element_count = 0
    parsed_formula = chemparse.parse_formula(formula)
    for k in parsed_formula.keys():
        sum_valence= sum_valence+parsed_formula[k]*valence_dict[k]
        element_count = element_count+parsed_formula[k]
        if valence_dict[k]>max_valence:
            max_valence=valence_dict[k]
    # if sum_valence%2 != 0:
    #     return False
    if sum_valence<2*max_valence:
        return False
    if sum_valence<2*(element_count-1):
        return False
    return True
def check_ratio(formula):
    """
    Checks the composition of chemical formula using ratio checks in the "7 golden rules" 
    ('Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry').

    Args:
        formula (str): The chemical formula to be checked.
    Returns:
        bool: True if the formula passes all checks, False otherwise.
        np.NAN: If the formula is invalid due to non-alphanumeric characters at the end.
    The function performs the following checks:
        - Checks the number of hydrogen and carbon atoms based on the accurate mass.
        - Checks the number of nitrogen and oxygen atoms.
        - Ensures the it is not a pure carbon/nitrogen loss (except N2)
        - Checks the hydrogen to carbon ratio.
        - Checks the fluorine to carbon ratio.
        - Checks the chlorine to carbon ratio.
        - Checks the bromine to carbon ratio.
        - Checks the nitrogen to carbon ratio.
        - Checks the oxygen to carbon ratio.
        - Checks the phosphorus to carbon ratio.
        - Checks the sulfur to carbon ratio.
        - Checks the silicon to carbon ratio.
    
    """

    if len(formula)==0:
        return False
    if formula[-1].isalnum() ==False:
        print(f'the formula passes {formula} is not right')
        return(np.NAN)
    parsed_formula = chemparse.parse_formula(formula)   
    accurate_mass = Formula(formula).isotope.mass
    if accurate_mass >0 and accurate_mass<500:
        if 'H' in parsed_formula.keys() and parsed_formula['H']>72:
            return False
        if 'C' in parsed_formula.keys() and parsed_formula['C']>39:
            return False
    elif accurate_mass>500 and accurate_mass<1000:
        if 'H' in parsed_formula.keys() and parsed_formula['H']>126:
            return False
        if 'C' in parsed_formula.keys() and parsed_formula['C']>78:
            return False

    if 'N' in parsed_formula.keys() and parsed_formula['N']>20:
        return False

    if 'O' in parsed_formula.keys() and parsed_formula['O']>27:
        return False
    if len(parsed_formula.keys())==1 and next(iter(parsed_formula.keys()))!=2:
        #check for pure carbon loss
        return(False)
    if len(parsed_formula.keys())==1 and next(iter(parsed_formula.keys()))=='N' and next(iter(parsed_formula.values()))=='N':
        #check for pure nitrogen loss (while not N2)
        return (False)
    if 'C' in parsed_formula.keys() and 'H' in parsed_formula.keys():
        if parsed_formula['H']/parsed_formula['C']>6 or parsed_formula['H']/parsed_formula['C']<0.1:
            # 7 golden rules: HC check
            return False
    if 'C' in parsed_formula.keys() and 'F' in parsed_formula.keys():
        if parsed_formula['F']/parsed_formula['C']>6:
            # 7 golden rules: CF check
            return False

    if 'C' in parsed_formula.keys() and 'Cl' in parsed_formula.keys():
        if parsed_formula['Cl']/parsed_formula['C']>2:
            # 7 golden rules: CCl check
            return False

    if 'C' in parsed_formula.keys() and 'Br' in parsed_formula.keys():
        if parsed_formula['Br']/parsed_formula['C']>2:
            # 7 golden rules: CF check
            return False

    if 'C' in parsed_formula.keys() and 'N' in parsed_formula.keys():
        if parsed_formula['N']/parsed_formula['C']>4:
            # 7 golden rules: CF check
            return False

    if 'C' in parsed_formula.keys() and 'O' in parsed_formula.keys():
        if parsed_formula['O']/parsed_formula['C']>3:
            # 7 golden rules: CF check
            return False

    if 'C' in parsed_formula.keys() and 'P' in parsed_formula.keys():
        if parsed_formula['P']/parsed_formula['C']>2:
            # 7 golden rules: CF check
            return False

    if 'C' in parsed_formula.keys() and 'S' in parsed_formula.keys():
        if parsed_formula['S']/parsed_formula['C']>3:
            # 7 golden rules: CF check
            return False

    if 'C' in parsed_formula.keys() and 'Si' in parsed_formula.keys():
        if parsed_formula['Si']/parsed_formula['C']>1:
            # 7 golden rules: CF check
            return False
    return(True)