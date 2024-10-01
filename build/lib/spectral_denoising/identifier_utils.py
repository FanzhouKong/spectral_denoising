import numpy as np 
from rdkit import Chem
import re
import requests
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from pubchempy import Compound, get_compounds
def get_classyfire(smiles, if_np=False):
    url = create_classyfire_url(smiles, if_np)
    r = requests.get(url)
    if r.ok:
        return r.json()
    else:
        return np.NAN
def everything_to_smiles(input):
    if input != input:# check for nan
        return np.NAN
    if is_smiles(input):
        smiles = input
    elif is_mol(input):
        smiles = Chem.MolToSmiles(input)

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
def everything_to_formula(input):
    if input != input:
        return np.NAN
    if is_smiles(input) == False:
        if is_formula(input):
            return input
    smiles = everything_to_smiles(input)
    mol = Chem.MolFromSmiles(smiles)
    formula_temp = CalcMolFormula(mol)
    # formula = standarize_formula(formula_temp)
    return(formula_temp)

def create_classyfire_url(smiles_string, if_np = True):
    if if_np:
        url_template = "https://npclassifier.gnps2.org/classify?smiles={}"
    else:
        url_template='https://structure.gnps2.org/classyfire?smiles={}'
    return url_template.format(smiles_string)

# Example usage

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
        return mol
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