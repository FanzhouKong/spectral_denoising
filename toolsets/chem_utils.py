from rdkit import Chem
import re
import chemparse
from molmass import Formula
from pubchempy import Compound, get_compounds
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Descriptors import ExactMolWt
import requests
import cirpy# reference plz see https://www.resources.aropha.com/blog/get-chemical-smiles-by-cas-or-name/
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdFMCS

import requests


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


def calculate_similarity(mol1, mol2):
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




def create_classyfire_url(smiles_string, if_np = True):
    if if_np:
        url_template = "https://npclassifier.gnps2.org/classify?smiles={}"
    else:
        url_template='https://structure.gnps2.org/classyfire?smiles={}'
    return url_template.format(smiles_string)

# Example usage
def get_classyfire(smiles, if_np=False):
    url = create_classyfire_url(smiles, if_np)
    r = requests.get(url)
    if r.ok:
        return r.json()
    else:
        return np.NAN
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
            uncharged_smiles = remove_acidic_hydrogen(components[0])
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
            uncharged_smiles = remove_acidic_hydrogen(uncharged)
        else:
            uncharged_smiles = Chem.MolToSmiles(uncharged)

        return( uncharged_smiles)
def remove_acidic_hydrogen(molecule, is_smiles = False):
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
def everything_to_formula(input):
    if input != input:
        return np.NAN
    smiles = everything_to_smiles(input)
    mol = Chem.MolFromSmiles(smiles)
    formula_temp = CalcMolFormula(mol)
    # formula = standarize_formula(formula_temp)
    return(formula_temp)
# def standarize_formula(formula):
#     if formula[-1] in ['+', '-']:
#         main_body = formula[0:-1]
#         polarity = formula[-1]
#     else:
#         main_body = formula
#         polarity = ''
#     formula_s = Formula(main_body).formula
#     formula_s = formula_s+polarity
#     return formula_s
def everything_to_smiles(input):
    if input != input:# check for nan
        return np.NAN
    if is_mol(input):
        smiles = Chem.MolToSmiles(input)
    elif is_smiles(input):
        smiles = input
    elif is_inchikey(input):
        smiles = inchikey_to_smiles(input)
    elif is_cas_number(input):
        smiles = cas_to_smiles(input)
        if smiles != smiles:
            smiles = name_to_smiles(input)
    else:
        smiles = name_to_smiles(input)
    return(smiles)
def everything_to_inchikey(input, first_block = True):
    smiles = np.NAN
    if input != input:# check for nan
        return np.NAN
    if is_inchikey(input):
        if first_block ==True:
            return input[0:14]
        else:
            return input
    elif is_mol(input):
        smiles = Chem.MolToSmiles(input)
    elif is_smiles(input):
        smiles = input
    elif is_cas_number(input):
        smiles = cas_to_smiles(input)
        if smiles != smiles:
            smiles = name_to_smiles(input)
    else:
        smiles = name_to_smiles(input)
    if smiles == smiles:
        mol = Chem.MolFromSmiles(smiles)
        inchikey = Chem.MolToInchiKey(mol)
        if first_block ==True:
            return inchikey[0:14]
        return inchikey
    else:
        return np.NAN
#below are individual building blocks to smiles
def smiles_to_inchikey(smiles):
    if isinstance(smiles, float):
        return np.NAN
    mol = Chem.MolFromSmiles(smiles)
    inchikey = Chem.MolToInchiKey(mol)
    return(inchikey[0:14])
def inchikey_to_smiles(inchikey):
    cc = get_compounds(inchikey, 'inchikey')
    if len(cc)>0:
        return (cc[0].isomeric_smiles)
    else:
        cc = get_compounds(inchikey[0:14], 'inchikey')
        if len(cc)>0:
            return (cc[0].isomeric_smiles)
        else:
            return (np.NAN)
def cas_to_smiles(cas):
    smile = cirpy.resolve(cas, 'smiles')
    if smile is None:
        smile = np.NAN
    return(smile)
def name_to_smiles(name):
    cc = get_compounds(name, 'name')
    if len(cc)>0:
        return (cc[0].isomeric_smiles)
    else:
        return (np.NAN)

def everything_to_image(molecule, savepath):
    from rdkit import Chem
    from rdkit.Chem import Draw
    if is_mol(molecule):
    # Create an RDKit molecule object
        mol = molecule

    elif is_smiles(molecule):
        # print('ttt')
        mol = Chem.MolFromSmiles(molecule)
    else:
        smiles = everything_to_smiles(molecule)
        mol = Chem.MolFromSmiles(smiles)
    # Generate the image of the molecule
    img = Draw.MolToImage(mol)
        # Save the image to a file
    img.save(savepath)
def break_adduct(adduct):
    if adduct != adduct:
        return np.NAN, np.NAN, np.NAN
    if '(i)' in adduct:
        return  np.NAN, np.NAN, np.NAN
    if adduct[0]=='[':
        adduct_part = re.split(r'\[|]', adduct)[1]
        charge = re.split(r'\[|]', adduct)[-1]
        ind = re.split(r'(\+|-)', adduct_part)
    else:
        adduct_part = adduct
        charge = re.split(r'\[|]', adduct)[-1]
        ind = re.split(r'(\+|-)', adduct_part)
    if charge[-1]=='+':
        polarity = 1
    elif charge[-1]=='-':
        polarity = -1
    else:
        return np.NAN, np.NAN, np.NAN
    if len(charge)>1:
        polarity = polarity*np.int32(charge[:-1])
    return(adduct_part, ind, polarity)
def parse_ind(s):
    #this is to parse individual instance of adduct
    index = 0
    coef = 1
    # Loop through each character in the string
    for char in s:
        # If the character is a digit, increment the index
        if char.isdigit():
            index += 1
        # If a non-digit character is found, break out of the loop
        else:
            break
    if index >0:
        coef = int(s[:index])
    # Return the part of the string without the leading integers
    return coef,s[index:]
def calculate_precursormz(formula, adduct, if_smiles = False):
    if adduct in ['Cat', 'CAT','[M+]', 'M+','[Cat]+']:
        adduct = '[M]+'
    if 'Hac' in adduct:
        adduct = adduct.replace("Hac", "C2H4O2")
    if 'FA' in adduct:
        adduct = adduct.replace("FA", "CH2O2")
    if 'DMSO' in adduct:
        adduct = adduct.replace("DMSO", "C2H6OS")
    if 'ACN' in adduct:
        adduct = adduct.replace("ACN", "C2H3N")
    if adduct[-1] in ['+','t']:
        adduct_polarity = 1
    elif adduct[-1]=='-':
        adduct_polarity = -1
    else:
        # print('cannot determine adduct polarity!')
        return (np.NAN)
    # if adduct[0]=='[':
    #     adduct_part = re.split(r'\[|]', adduct)[1]
    #     ind = re.split(r'(\+|-)', adduct_part)
    # else:
    #     adduct_part = adduct
    #     ind = re.split(r'(\+|-)', adduct_part)
    adduct_part, ind, charge = break_adduct(adduct)
    coef,a  = parse_ind(ind[0])

    if if_smiles == True:
        mol = Chem.MolFromSmiles(formula)
        # max_c_length = find_longest_element_chain(mol, 'C')
        formula = CalcMolFormula(mol)
    if formula[-1] in ['+','-'] and adduct_part not in ['M']:

        return 0
    elif formula[-1] not in ['+', '-'] and adduct_part in ['M']:
        return 0
    elif formula[-1] in ['+', '-'] and adduct_part in ['M']:
        if formula[-1]==adduct[-1]:
            accurate_mass = Formula(formula[:-1]).isotope.mass
            accurate_mass = accurate_mass-adduct_polarity*0.00054857990924
            return accurate_mass
        else:
            return 0



    accurate_mass = Formula(formula).isotope.mass

    accurate_mass = accurate_mass*coef


    for i in range(1, len(ind)):
        if (ind[i] not in ['+', '-']) and len(ind[i])>0:
            coef, a = parse_ind(ind[i])
            if ind[i-1]=='+':
                polarity = 1
            elif ind[i-1]=='-':
                polarity = -1
            else:
                polarity = 0

            accurate_mass = accurate_mass+polarity*coef*Formula(a).isotope.mass
    accurate_mass = accurate_mass-adduct_polarity*0.00054857990924
    accurate_mass = accurate_mass/abs(charge)
    return accurate_mass


#below are is_ section
def is_inchikey(string):
    # Define the regex pattern for InChIKeys
    pattern = r'^[A-Z]{14}-[A-Z]{10}-[A-Z0-9]$'

    # Use re.match to check if the pattern matches the entire string
    if re.match(pattern, string):
        return True
    else:
        return False

def is_mol(obj):
    return isinstance(obj, Chem.rdchem.Mol)
def is_smiles(smiles_string):
    # Attempt to create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_string)

    # If the molecule object is created successfully, the SMILES string is valid
    if mol is not None:
        return True
    else:
        # If the molecule object is None, the SMILES string is invalid
        return False
def is_cas_number(string):
    # Regex pattern for CAS numbers: one or more digits, followed by a hyphen, followed by two or more digits,
    # followed by a hyphen, and ending with a single digit
    pattern = r'^\d+-\d{2,}-\d$'

    # Check if the given string matches the pattern
    if re.match(pattern, string):
        return True
    else:
        return False
def everything_to_mw(mol):
    if is_mol(mol)==False:
        smiles = everything_to_smiles(mol)
        mol = Chem.MolFromSmiles(smiles)
    return(ExactMolWt(mol))
def is_formula(s):
    # Regular expression to match chemical formulas
    # Starts with an uppercase letter, optionally followed by a lowercase letter (for two-letter elements)
    # Optionally followed by a number (for the count of atoms)
    # This pattern repeats throughout the string
    pattern = r'^([A-Z][a-z]?\d*)+$'

    # Match the entire string against the pattern
    match = re.fullmatch(pattern, s)

    # If there's a match, the string is a valid chemical formula
    return bool(match)