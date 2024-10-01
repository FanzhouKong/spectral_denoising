
__all__ = ['get_bond_similarity','__remove_bonds_in_smarts', 'desalter','__remove_acidic_hydrogen','calculate_precursormz']

from rdkit import Chem
import re
import chemparse
from molmass import Formula

from rdkit import Chem

from rdkit.Chem.MolStandardize import rdMolStandardize

import requests
import cirpy# reference plz see https://www.resources.aropha.com/blog/get-chemical-smiles-by-cas-or-name/
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdFMCS

import requests

from .identifier_utils import *


def parse_formula(formula):
    dict = chemparse.parse_formula(formula)
    lst = ([[x,y] for x, y in zip([*dict],[*dict.values()]) ])
    return(lst)

def transpose_formula(lst):
    
    return([list(x) for x in zip(*lst)])
def get_bond_similarity(mol1, mol2):
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
        # print(Chem.MolToSmiles(mol1))
        # print(Chem.MolToSmiles(mol2))
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
        # print(common_s)

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
def __remove_acidic_hydrogen(molecule, is_smiles = False):
    # Convert the SMILES string to a RDKit molecule object
    if is_smiles==True:
        smiles = molecule
        molecule = Chem.MolFromSmiles(molecule)
    else:
        smiles = Chem.MolToSmiles(molecule)
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


def calculate_precursormz(adduct_string,mol = None, testing = False):
    if testing == True:
        molecule_mass = 853.33089
    else:
        molecule_mass = Formula(everything_to_formula(mol)).isotope.mass
    
    adduct_string = replace_adduct_string(adduct_string)
    electron_mass = 0.00054858026
    m_coef = determine_parent_coefs(adduct_string)
    charge = determine_adduct_charge(adduct_string)
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
    m_pattern = r'(\d*)M'
    m_match = re.search(m_pattern, adduct_string)
    if m_match:
        # If the match is empty (no number before 'M'), default the coefficient to 1
        coefficient = m_match.group(1)
        coefficient = int(coefficient) if coefficient else 1
        return coefficient
    else:
        print(f'the correct adduct form cannot be determined from {adduct_string}')
        return np.NAN
def determine_adduct_charge(adduct_string):
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
        charge = np.NAN
        print(f'the correct adduct form cannot be determined from {adduct_string}')
    return charge
def parse_adduct(adduct_string):
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
    if adduct_string in ['Cat', 'CAT','[M+]', 'M+','[Cat]+']:
        adduct_string = '[M]+'
    if adduct_string in ['Cat-', 'CAT-','[M]-', 'M-','[Cat]-']:
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