import chemparse
import numpy as np
from molmass import Formula
from .constant import valence_dict
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
    # if sum_valence%2 != 0:
    #     return False
    if sum_valence<2*max_valence:
        return False
    if sum_valence<2*(element_count-1):
        return False
    return True
def check_ratio(formula):
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