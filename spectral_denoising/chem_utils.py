

from rdkit import Chem
import re
import chemparse
from molmass import Formula

from rdkit import Chem

from rdkit.Chem.MolStandardize import rdMolStandardize

import requests
import cirpy
# reference plz see https://www.resources.aropha.com/blog/get-chemical-smiles-by-cas-or-name/
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdFMCS

import requests

from .identifier_utils import *

def get_bond_similarity(mol1, mol2):
    """
    Calculate the bond similarity between two molecules.
    Detailed algorithm can be found in the following paper: Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification (Yuanyue Li et al., 2021)

    Args:
        mol1: The first molecule, which can be in various formats (e.g., SMILES string, RDKit molecule object).
        mol2: The second molecule, which can be in various formats (e.g., SMILES string, RDKit molecule object).
    Returns:
        dict: A dictionary containing the following keys:
            - "mol1_bond_number": The number of bonds in the first molecule.
            - "mol2_bond_number": The number of bonds in the second molecule.
            - "common_bond_number": The number of common bonds between the two molecules.
            - "bond_difference": The average bond difference between the two molecules.
            - "bond_similarity": The bond similarity score between the two molecules, calculated as 
              (2 * common_bond_number) / (mol1_bond_number + mol2_bond_number).
            - "minimal_diff": The minimal difference in bond numbers after removing common bonds.
    """
    
    if is_mol(mol1)==False:
        smile1 = everything_to_smiles(mol1)
        mol1 = Chem.MolFromSmiles(smile1)
    if is_mol(mol2)==False:
        smile2 = everything_to_smiles(mol2)
        mol2 = Chem.MolFromSmiles(smile2)
    mol1 = Chem.rdchem.RWMol(mol1)
    mol2 = Chem.rdchem.RWMol(mol2)
    bond_number_common = 0
    bond_number_mol1 = len(mol1.GetBonds())
    bond_number_mol2 = len(mol2.GetBonds())
    while True:
        res = rdFMCS.FindMCS(
            [mol1, mol2],
            timeout=20,
            threshold=1,
            # ringMatchesRingOnly=True,
            # completeRingsOnly=True,
            # atomCompare=rdFMCS.AtomCompare.CompareElements,
            # bondCompare=rdFMCS.BondCompare.CompareOrderExact,
            # ringCompare=rdFMCS.RingCompare.StrictRingFusion,
            # maximizeBonds=True,
            # matchValences=True,
        )
        if res.numBonds == 0:
            break

        common_s = res.smartsString
        mol1, _ = __remove_bonds_in_smarts(mol1, common_s)
        mol2, _ = __remove_bonds_in_smarts(mol2, common_s)
        bond_number_common += res.numBonds
        # print(bond_number_common)
        minimal_diff = np.min([bond_number_mol1-bond_number_common, bond_number_mol2-bond_number_common])
    return {
        "mol1_bond_number": bond_number_mol1,
        "mol2_bond_number": bond_number_mol2,
        "common_bond_number": bond_number_common,
        "bond_difference": (bond_number_mol1 + bond_number_mol2) / 2 - bond_number_common,
        "bond_similarity": (2 * bond_number_common) / (bond_number_mol1 + bond_number_mol2),
        'minimal_diff':minimal_diff
    }

def __remove_bonds_in_smarts(mol, smarts):
    """
    Helper function for get_bond_similarity
    Remove bonds in a molecule based on a SMARTS pattern.
    This function takes a molecule and a SMARTS pattern, identifies the substructure
    in the molecule that matches the SMARTS pattern, and removes the bonds between
    atoms in the matched substructure.
    Args:
        mol (rdkit.Chem.Mol): The molecule from which bonds will be removed.
        smarts (str): The SMARTS pattern used to identify the substructure.
    Returns:
        tuple: A tuple containing the modified molecule (rdkit.Chem.Mol) and the number
               of bonds removed (int).
    """

    removed_bond_number = 0
    pattern = Chem.MolFromSmarts(smarts)
    Chem.GetSymmSSSR(pattern)
    Chem.GetSymmSSSR(mol)
    sub_atom = mol.GetSubstructMatch(pattern)
    # print(sub_atom)
    for i in sub_atom:
        all_bonds = mol.GetAtomWithIdx(i).GetBonds()
        for bond in all_bonds:
            atom_1, atom_2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            mol.RemoveBond(atom_1, atom_2)
            if atom_1 in sub_atom and atom_2 in sub_atom:
                removed_bond_number += 1
    return mol, removed_bond_number

def desalter(input):

    """
    Processes the input molecule to remove salts and return an uncharged SMILES string.

    Args:
        input (str or RDKit Mol): The input molecule, which can be a SMILES string or an RDKit Mol object.
    Returns:
        uncharged_smiles (str): The uncharged SMILES string of the largest component of the input molecule.
    Notes:
        - If the input is not a valid molecule, the function will attempt to convert it to a SMILES string.
        - If the input molecule contains multiple components, the largest component will be processed.
        - If the largest component has a formal charge of +1, acidic hydrogens will be removed.
        - If the input is NaN, the function will return np.NAN.
    """

    if input != input:
        return np.NAN
    if is_mol(input) == False:
        smile = everything_to_smiles(input)
        mol = Chem.MolFromSmiles(smile)
    else:
        mol = input
    components = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

    if len(components)==1:
        if Chem.GetFormalCharge(components[0])==1:
            uncharged_smiles = __remove_acidic_hydrogen(components[0])
        else:
            uncharged_smiles = Chem.MolToSmiles(components[0])
        return uncharged_smiles
    else:
        n_atoms = np.zeros(len(components))
        counter = 0
        for component in components:
            n_atoms[counter]=component.GetNumAtoms()
            counter = counter+1
        idx = np.argmax(n_atoms)
        # print(idx)
        charged = components[np.argmax(n_atoms)]
        un = rdMolStandardize.Uncharger()
        uncharged = un.uncharge(charged)
        # return(uncharged)
        if Chem.GetFormalCharge(uncharged)==1:
            uncharged_smiles = __remove_acidic_hydrogen(uncharged)
        else:
            uncharged_smiles = Chem.MolToSmiles(uncharged)

        return( uncharged_smiles)
def __remove_acidic_hydrogen(mol):
    """
    Helper function of desalter
    Remove the acidic hydrogen from carboxylic acid groups in a molecule.
    This function identifies carboxylic acid groups in a given molecule and removes the hydrogen atom from the hydroxyl group, 
    setting the formal charge of the oxygen atom to -1. The input molecule can be provided either as an RDKit molecule object 
    or as a SMILES string.
    Args:
        molecule (Union[rdkit.Chem.rdchem.Mol, str]): The input molecule, either as an RDKit molecule object or a SMILES string.
        is_smiles (bool, optional): A flag indicating whether the input molecule is a SMILES string. Defaults to False.
    Returns:
        str: The SMILES string of the modified molecule with the acidic hydrogen removed. If no carboxylic acid group is found, 
        the original SMILES string is returned.
    """
    # Convert the SMILES string to a RDKit molecule object
    if is_mol(mol)==False:
        smiles = everything_to_smiles(mol)
        molecule = Chem.MolFromSmiles(smiles)
    else:
        molecule = mol
    # Define the SMARTS pattern for carboxylic acids (includes the hydrogen in the hydroxyl group)
    carboxylic_acid_smarts = 'C(=O)[OH]'
    # Create a query molecule from the SMARTS pattern
    query = Chem.MolFromSmarts(carboxylic_acid_smarts)
    # Find substructures that match the query (carboxylic acids)
    matches = molecule.GetSubstructMatches(query)
    if not matches:
        # print("No carboxylic acid group found.")
        return smiles  # Return the original SMILES if no carboxylic acid group is found
    editable_mol = Chem.RWMol(molecule)
    # Assuming only one carboxylic group needs to be modified,
    # and focusing on the first match
    for match in matches:
        # The oxygen atom in the OH group is the last second in the matched pattern
        oxygen_idx = match[-1]
        # Get the oxygen atom
        oxygen_atom = editable_mol.GetAtomWithIdx(oxygen_idx)
        # Set the formal charge of the oxygen atom to -1
        oxygen_atom.SetFormalCharge(-1)

        # Set the implicit hydrogen count of the oxygen atom to 0
        # Assuming there's only one hydrogen bonded which we want to remove
        oxygen_atom.SetNumExplicitHs(0)

        # Break after the first modification, assuming only one modification is needed
        break

    # Convert back to a molecule
    modified_mol = editable_mol.GetMol()

    # Convert the modified molecule back to SMILES without sanitization
    modified_smiles = Chem.MolToSmiles(modified_mol)

    return modified_smiles
def parse_formula(formula):

    """
    Parses a chemical formula into its constituent elements and their quantities.

    Args:
        formula (str): A string representing the chemical formula (e.g., "H2O", "C6H12O6").
    Returns:
         list: A list of lists, where each inner list contains an element and its quantity 
         (e.g., [['H', 2], ['O', 1]] for "H2O").
    """

    dict = chemparse.parse_formula(formula)
    lst = ([[x,y] for x, y in zip([*dict],[*dict.values()]) ])
    return(lst)

def transpose_formula(lst):

    """
    Transpose a parsed formula in a nested list format from [[element, quantity], ...] to [[element, ...], [quantity, ...]].

    Args:
        lst (list): A list of lists where each sublist represents a row of the matrix.
    Returns:
        list: A transposed version of the input list of lists, where rows are converted to columns and vice versa.
    """
    
    return([list(x) for x in zip(*lst)])
def calculate_precursormz(adduct_string,mol = None, testing = False):
    
    """
    Calculate the precursor m/z (mass-to-charge ratio) for a given molecule and adduct string.
    Very robust function, handles a wide variety of adducts strings.

    Args:
        adduct_string (str): The adduct string representing the ion type and charge state.
        mol (str, optional): The molecular formula of the compound. Defaults to None.
        testing (bool, optional): If True, use a predefined molecule mass for testing purposes. Defaults to False.
        float: The calculated precursor m/z value.
    Returns:
        precursor_mz (float): The calculated precursor m/z value.
    Raises:
        Warning: If an unrecognized adduct is encountered in the adduct string, it will be ignored and a warning will be printed.
    Notes:
        - The function uses the `Formula` class to calculate the mass of the molecule and ions.
        - The `replace_adduct_string`, `determine_parent_coefs`, `determine_adduct_charge`, and `parse_adduct` functions are assumed to be defined elsewhere in the codebase.
        - The electron mass is considered in the calculation to adjust for the loss/gain of electrons.
    """
    
    if testing == True:
        molecule_mass = 853.33089
    else:
        molecule_mass = Formula(everything_to_formula(mol)).isotope.mass

    adduct_string = adduct_string.strip()
    adduct_string = replace_adduct_string(adduct_string)
    if 'i' in adduct_string or 'M' not in adduct_string:
        return np.nan
    electron_mass = 0.00054858026
    m_coef = determine_parent_coefs(adduct_string)
    charge = determine_adduct_charge(adduct_string)
    if m_coef != m_coef or charge != charge:
        return np.nan
    # Check for charge in the adduct string
    
        

    # Find all matches of adduct parts
    adduct = parse_adduct(adduct_string)
    mass_change = 0.0
    for adduct in adduct:
        sign, count, ion_type = adduct
        # Determine the sign (+ or -)
        sign_multiplier = 1 if sign == '+' else -1
        try:
            ion_mass = Formula(ion_type).isotope.mass
        except:
            print(f'Warning: Unrecognized adduct {ion_type} in {adduct_string} is ignored.')
            continue
        mass_change += sign_multiplier * count * ion_mass
        # Calculate the mass change if ion_type is recognized
    mass_change = mass_change-charge*electron_mass
    # print(mass_change)

    # Calculate the precursor m/z, considering the charge state, adjust for loss/gain of electrons
    precursor_mz = molecule_mass*m_coef+mass_change
    precursor_mz = precursor_mz/ abs(charge)
    return precursor_mz
def determine_parent_coefs(adduct_string):
    """
    Determine the coefficient of of the adducts by ignoring the parent ion M.

    Args:
        adduct_string (str): The adduct string from which to determine the parent molecule coefficient.
    Returns:
        coefficient (int): The coefficient of the adduct (e.g. M+H, coef =1, M+2H+, coef = 2). If the adduct string does not match the expected pattern, the function prints an error message and returns numpy's NaN.
    """

    m_pattern = r'(\d*)M'
    m_match = re.search(m_pattern, adduct_string)
    if m_match:
        # If the match is empty (no number before 'M'), default the coefficient to 1
        coefficient = m_match.group(1)
        coefficient = int(coefficient) if coefficient else 1
        return coefficient
    else:
        # print(f'the correct adduct form cannot be determined from {adduct_string}')
        return np.nan
def determine_adduct_charge(adduct_string):
    """
    Determine the charge of an adduct based on its string representation.
    This function processes an adduct string to determine its charge. The adduct string
    is first standardized using the `replace_adduct_string` function. The charge is then
    determined based on the ending character(s) of the string.

    Args:
        adduct_string (str): The string representation of the adduct.
    Returns:
        charge (int): The charge of the adduct. Returns a positive integer for positive charges, 
        a negative integer for negative charges, and NaN if the charge cannot be determined.
    Notes:
        - If the adduct string ends with '+', the function checks if the preceding character is a digit
          and if the charge is enclosed in brackets. If so, it extracts the charge; otherwise, it assumes
          a charge of +1.
        - If the adduct string ends with '-', the function performs a similar check for negative charges.
        - If the adduct string does not end with '+' or '-', the function returns NaN and prints a message
          indicating that the charge could not be determined.
    """

    adduct_string = replace_adduct_string(adduct_string)
    if adduct_string.endswith('+'):
        if adduct_string[-2].isdigit() and adduct_string[-3] == ']': # Check if the charge is in brackets
            charge = int(adduct_string[-2])
        else:
            charge = 1
    elif adduct_string.endswith('-'):
        if adduct_string[-2].isdigit() and adduct_string[-3] == ']': # Check if the charge is in brackets:
            charge = -int(adduct_string[-2])
        else:
            charge = -1
    else:
        charge = np.nan
        # print(f'the correct adduct form cannot be determined from {adduct_string}')
    return charge
def parse_adduct(adduct_string):
    """
    Parses an adduct string into its components.
    This function takes an adduct string and breaks it down into its constituent parts,
    including the sign, count, and ion type. The adduct string is first processed to 
    replace certain patterns, and then a regular expression is used to capture the 
    different parts of the adduct.

    Args:
        adduct_string (str): The adduct string to be parsed.
    Returns:
        list: A list of lists, where each sublist contains the sign (str), count (int), and ion type (str) of each part of the adduct.
    """

    adduct_string = replace_adduct_string(adduct_string)
    # Regular expression to capture multiple parts of the adduct and the charge
    pattern = r'([+-])(\d*)([A-Za-z0-9]+)'
    matches = re.findall(pattern, adduct_string)
    matches=map(list, matches)
    parsed_adduct = []
    for match in matches:
            sign, count, ion_type = match
            
            # If no count is given, assume it's 1
            count = int(count) if count else 1
            parsed_adduct.append([sign, count, ion_type])
            # Determine the sign (+ or -)
    return parsed_adduct
def replace_adduct_string(adduct_string):
    """
    Replaces specific adduct strings with their standardized or chemical formula equivalents.
    This function takes an adduct string and replaces it with a standardized version or its 
    corresponding chemical formula.

    Args:
        adduct_string (str): The adduct string to be replaced.
    Returns:
        adduct_string (str): The replaced adduct string.
    """
    adduct_string=adduct_string.replace('Cat', 'M')
    adduct_string=adduct_string.replace('CAT', 'M')
    if adduct_string in ['Cat', 'CAT','[M+]', 'M+','[Cat]+']:
        adduct_string = '[M]+'
    if adduct_string in ['[M]-', 'M-']:
        adduct_string = '[M]-'
    if 'Hac' in adduct_string:
        adduct_string = adduct_string.replace("Hac", "C2H4O2")
    if 'FA' in adduct_string:
        adduct_string = adduct_string.replace("FA", "CH2O2")
    if 'DMSO' in adduct_string:
        adduct_string = adduct_string.replace("DMSO", "C2H6OS")
    if 'ACN' in adduct_string:
        adduct_string = adduct_string.replace("ACN", "C2H3N")
    return(adduct_string)