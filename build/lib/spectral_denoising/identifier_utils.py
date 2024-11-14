import numpy as np 
from rdkit import Chem
import re
import cirpy
import requests
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from pubchempy import Compound, get_compounds
from molmass import Formula
def get_classyfire(smiles, if_np=False):
    """
    Retrieves the ClassyFire classification for a given SMILES string.
    
    Args:
        smiles (str): The SMILES string of the molecule to classify.
        if_np (bool, optional): A flag indicating whether the molecule is a natural product. Defaults to False.
    Returns:
        dict: The JSON response from the ClassyFire API if the request is successful, 
                       otherwise numpy.NAN.
    """

    url = create_classyfire_url(smiles, if_np)
    r = requests.get(url)
    if r.ok:
        return r.json()
    else:
        return np.nan
def everything_to_smiles(input):
    """
    Convert various chemical identifier formats to a SMILES string.
    This function takes an input which can be in different chemical identifier formats 
    (SMILES, Mol, InChIKey, CAS number, or chemical name) and converts it to a SMILES string.
    
    Args:
        input (str or RDKit Mol): The chemical identifier to be converted. It can be a SMILES string, an RDKit Mol object, an InChIKey, a CAS number, or a chemical name.
    Returns:
        smiles (str): The corresponding SMILES string if the conversion is successful. Returns NaN if the input is NaN.
    """

    if input != input:# check for nan
        return np.nan
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
    """
    Converts various chemical identifiers to an InChIKey or its first block.
    This function takes an input which can be an InChIKey, a molecule object, a SMILES string, 
    a CAS number, or a chemical name, and converts it to an InChIKey. If the input is already 
    an InChIKey, it can return either the full InChIKey or just the first block of it based 
    on the `first_block` parameter.
    
    Args:
        input (str or RDKit Mol): The chemical identifier to be converted. It can be an InChIKey, a molecule object, a SMILES string, a CAS number, or a chemical name.
        first_block (bool, optional): If True, returns only the first block of the InChIKey. Defaults to True.
    Returns:
        inchikey (str): The InChIKey or its first block if `first_block` is True. Returns NaN if the 
        input is invalid or cannot be converted.
    """

    smiles = np.nan
    if input != input:# check for nan
        return np.nan
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
        return np.nan
def everything_to_formula(input):
    """
    Converts various chemical input to a molecular formula.
    This function takes an input which can be in different chemical formats 
    (e.g., SMILES, molecular formula) and converts it to a standardized 
    molecular formula. If the input is already a molecular formula, it is 
    returned as is. If the input is a SMILES string, it is first converted 
    to a molecular object and then to a molecular formula.

    Args:
        input (str): The chemical input which can be a SMILES string, 
        molecular formula, or other recognizable chemical format.
    Returns:
        formula (str): The molecular formula of the input chemical. If the input is 
             invalid or cannot be converted, returns NaN.
    """

    if input != input:
        return np.nan
    if is_smiles(input) == False:
        if is_formula(input):
            return input
    smiles = everything_to_smiles(input)
    mol = Chem.MolFromSmiles(smiles)
    formula = CalcMolFormula(mol)
    # formula = standarize_formula(formula_temp)
    return(formula)

def create_classyfire_url(smiles_string, if_np = True):
    """
    Generates a URL for ClassyFire or NPClassifier based on the provided SMILES string. Just a helper function
    """
    if if_np:
        url_template = "https://npclassifier.gnps2.org/classify?smiles={}"
    else:
        url_template='https://structure.gnps2.org/classyfire?smiles={}'
    return url_template.format(smiles_string)


def smiles_to_inchikey(smiles):
    """
    helper function

    Args:
        smiles (str): A SMILES string representing a molecule.
    Returns:
        inchikey (str): The InChIKey of the molecule, first block only.
    """
    if isinstance(smiles, float):
        return np.nan
    mol = Chem.MolFromSmiles(smiles)
    inchikey = Chem.MolToInchiKey(mol)
    return(inchikey[0:14])
def inchikey_to_smiles(inchikey):
    """
    helper function, but uses pubchem database

    Args:
        inchikey (str): The inchikey of the molecule to look up.
    Returns:
        str: The fetched isomeric SMILES code.
    """
    cc = get_compounds(inchikey, 'inchikey')
    if len(cc)>0:
        return (cc[0].isomeric_smiles)
    else:
        cc = get_compounds(inchikey[0:14], 'inchikey')
        if len(cc)>0:
            return (cc[0].isomeric_smiles)
        else:
            return (np.nan)
def cas_to_smiles(cas):
    """
    Convert a CAS (Chemical Abstracts Service) number to a SMILES (Simplified Molecular Input Line Entry System) string.

    Args:
        cas (str): The CAS number of the chemical compound.
    Returns:
        str: The SMILES string of the chemical compound if found, otherwise NaN.
    """

    smile = cirpy.resolve(cas, 'smiles')
    if smile is None:
        smile = np.nan
    return(smile)
def name_to_smiles(name):
    """
    Convert a chemical name to its corresponding SMILES (Simplified Molecular Input Line Entry System) representation, with Pubchem as backend.

    Args:
        name (str): The chemical name to be converted.
    Returns:
        str: The SMILES representation of the chemical if found, otherwise numpy.nan.
    """

    cc = get_compounds(name, 'name')
    if len(cc)>0:
        return (cc[0].isomeric_smiles)
    else:
        return (np.nan)

def everything_to_image(molecule, savepath):
    """
    Converts a molecular representation to an image and saves it to the specified path.

    Args:
        molecule (str or RDKit Mol object): The molecular representation, which can be a SMILES string, 
        an RDKit Mol object, or any other format that can be converted to a SMILES string.
        savepath (str): The file path where the generated image will be saved.
    Returns:
        None
    """

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
    """
    Check if a given string is a valid InChIKey using regex.
    An InChIKey is a 27-character string divided into three blocks by hyphens:
    - The first block contains 14 uppercase letters.
    - The second block contains 10 uppercase letters.
    - The third block contains a single uppercase letter or digit.

    Args:
        string (str): The string to be checked.
    Returns:
        bool: True if the string is a valid InChIKey, False otherwise.
    """

    # Define the regex pattern for InChIKeys
    pattern = r'^[A-Z]{14}-[A-Z]{10}-[A-Z0-9]$'
    # Use re.match to check if the pattern matches the entire string
    if re.match(pattern, string):
        return True
    else:
        return False

def is_mol(mol):
    """
    Check if the given object is an instance of Chem.rdchem.Mol.

    Args:
        mol: The object to check.
    Returns:
        bool: True if the object is an instance of Chem.rdchem.Mol, False otherwise.
    """

    return isinstance(mol, Chem.rdchem.Mol)
def is_smiles(smiles_string):
    """
    Check if a given string is a valid SMILES (Simplified Molecular Input Line Entry System) representation.

    Args:
        smiles_string (str): The SMILES string to be validated.
    Returns:
        bool: True if the SMILES string is valid, False otherwise.
    Example:
        >>> is_smiles("CCO")
        True
        >>> is_smiles("invalid_smiles")
        False
    """

    # Attempt to create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_string)
    # If the molecule object is created successfully, the SMILES string is valid
    if mol is not None:
        return True
    else:
        # If the molecule object is None, the SMILES string is invalid
        return False
def is_cas_number(string):
    """
    Check if a given string is a valid CAS (Chemical Abstracts Service) number.
    A CAS number is a unique numerical identifier assigned to every chemical substance
    described in the open scientific literature. It is formatted as one or more digits,
    followed by a hyphen, followed by two or more digits, followed by another hyphen,
    and ending with a single digit.

    Args:
        string (str): The string to be checked.
    Returns:
        bool: True if the string is a valid CAS number, False otherwise.
    """

    # Regex pattern for CAS numbers: one or more digits, followed by a hyphen, followed by two or more digits,
    # followed by a hyphen, and ending with a single digit
    pattern = r'^\d+-\d{2,}-\d$'

    # Check if the given string matches the pattern
    if re.match(pattern, string):
        return True
    else:
        return False
def everything_to_mw(mol):
    """
    Converts a given molecule representation to its molecular weight (MW).
    This function first checks if the input is a valid molecule object. If not, it attempts to convert the input to a SMILES string and then to a molecule object. Finally, it calculates and returns the exact molecular weight of the molecule.
    
    Args:
        mol: The input molecule representation. This can be a molecule object or another representation that can be converted to a SMILES string.
    Returns:
        float: The exact molecular weight of the molecule.
    Raises:
        ValueError: If the input cannot be converted to a valid molecule object.
    """
    
    if is_mol(mol)==False:
        smiles = everything_to_smiles(mol)
        
        mol = Chem.MolFromSmiles(smiles)
        return mol
    return(ExactMolWt(mol))
def is_formula(s):
    """
    Check if a given string is a valid chemical formula.
    A valid chemical formula starts with an uppercase letter, optionally followed by a lowercase letter 
    (for two-letter elements), and optionally followed by a number (for the count of atoms). This pattern 
    repeats throughout the string.

    Args:
        s (str): The string to be checked.
    Returns:
        bool: True if the string is a valid chemical formula, False otherwise.
    """

    # Regular expression to match chemical formulas
    # Starts with an uppercase letter, optionally followed by a lowercase letter (for two-letter elements)
    # Optionally followed by a number (for the count of atoms)
    # This pattern repeats throughout the string
    try:
        Formula(s).formula
        return True
    except:
        return False
    # pattern = r'^([A-Z][a-z]?\d*)+$'
    # if s.endswith('+') or s.endswith('-'):
    #     s = s[:-1]
    # # Match the entire string against the pattern
    # match = re.fullmatch(pattern, s)

    # # If there's a match, the string is a valid chemical formula
    # return bool(match)