from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import Draw
from spectral_denoising.spectral_denoising import has_benzene
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import re
import chemparse
from rdkit.Chem.Descriptors import ExactMolWt

from spectral_denoising.chem_utils import parse_formula, transpose_formula
def verify_status():
    print('yes spectral denoising is correctly installed')
def get_atom_index(mol, atom):
    lst = []
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetSymbol() == atom:
            lst.append(i)
    return lst

def get_rarest_heavy_atom(formula):
    """
    Get the rarest heavy atom in the dataset
    """
    parsed_formula = parse_formula(formula)
    atom, count = transpose_formula(parsed_formula)
    order = np.argsort(count)
    atom = np.array(atom)[order]
    count = np.array(count)[order]
    # print (atom, count)
    if atom[0] == 'H':
        return atom[1]
    else:
        return atom[0]
def enumaerate_candidate_list(frags, all_subformulas, all_sub_masses, mass_tolerance = 0.005):
    candidate_list = []

    for frag in frags:
        idx_left, idx_right = all_sub_masses.searchsorted([frag-proton_mass-mass_tolerance, frag-proton_mass+mass_tolerance])
        # print(idx_left, idx_right)
        candidate_list.append(all_subformulas[idx_left:idx_right])
    return candidate_list
def check_fragment_parallel(candidates, mol):
    # Use a process pool for parallelism
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Submit each candidate along with its index to preserve order
        future_to_index = {executor.submit(check_subformulas_topo, candidate, mol): idx for idx, candidate in enumerate(candidates)}
        
        # Collect results in a way that maintains the original order
        results = [None] * len(candidates)
        for future in concurrent.futures.as_completed(future_to_index):
            index = future_to_index[future]
            results[index] = future.result()
    
    return results
def check_fragment( all_candidates, mol):
    result = []
    for c in tqdm(all_candidates):
        result.append(check_subformulas_topo(c, mol))
    return result
def check_subformulas_topo(subformulas, mol):
    for subformula in subformulas:
        path_length = num_non_h_atoms(subformula)-1
        if check_topo(mol, subformula, int(path_length)):
            return True
    return False
def compare_formula(formula1: str, formula2: str) -> bool:
    modified_formula1 = re.sub(r'H\d*', '', formula1)
    modified_formula2 = re.sub(r'H\d*', '', formula2)
    return chemparse.parse_formula(modified_formula1) == chemparse.parse_formula(modified_formula2)
def check_topo(mol, formula, path_length: int):
    
    all_subgraphs = Chem.rdmolops.FindAllSubgraphsOfLengthN(mol, path_length)
    for subgraph in all_subgraphs:
        submol = Chem.PathToSubmol(mol, subgraph)
        sub_formula = CalcMolFormula(submol)
        if compare_formula(sub_formula, formula):
            return True
    return(False)
def get_subgraph_length(subformula):
    return int(num_non_h_atoms(subformula)-1)
def num_non_h_atoms(formula):
    formula_dict = chemparse.parse_formula(formula)
    sum_values = sum(value for key, value in formula_dict.items() if key != 'H')
    return int(sum_values)
# def break_bond_and_get_fragments(mol):
#     """Break each bond of the molecule and return the resulting fragments."""
#     Chem.Kekulize(mol,clearAromaticFlags=True)
#     bonds = mol.GetBonds()
#     broken_fragments = []
#     bond_idx = []
#     for bond in bonds:
#         idx_1 = bond.GetBeginAtomIdx()
#         idx_2 = bond.GetEndAtomIdx()

#         mol_copy = Chem.RWMol(mol)

#         mol_copy.RemoveBond(idx_1, idx_2)
#         fragments = rdmolops.GetMolFrags(mol_copy, asMols=True, sanitizeFrags=False)
#         broken_fragments.append(fragments)
#         bond_idx.append(bond.GetIdx())

#     return broken_fragments, bond_idx
# def break_ring_and_get_fragments(mol, ring_info):
#     pass
# def check_membership(bond_idx, bond_rings):
#     membership = [x for x in range(len(bond_rings)) if 42 in bond_rings[x]]
#     return membership
# def recursive_bond_break(mol):
#     Chem.Kekulize(mol,clearAromaticFlags=True)
#     """Recursively break the bonds of the molecule until each fragment pair has exactly two pieces."""
#     fragments = break_bond_and_get_fragments(mol)
    
#     result_fragments = []

#     for fragment_pair in fragments:
#         if len(fragment_pair) == 2:
#             # If already broken into two parts, add it to the result
#             result_fragments.append(fragment_pair)
#         elif len(fragment_pair) == 1:
#             # If only one fragment, keep breaking it recursively
#             sub_fragments = recursive_bond_break(fragment_pair[0])
#             result_fragments.extend(sub_fragments)

#     return result_fragments

# def parallel_recursive_bond_break(mol):
#     print('cleared aromatic flags')
#     Chem.Kekulize(mol,clearAromaticFlags=True)
#     """Parallelizes the recursive bond breaking process."""
#     fragments = break_bond_and_get_fragments(mol)
#     result_fragments = []
    
#     with ProcessPoolExecutor(6) as executor:
#         futures = []
        
#         for fragment_pair in fragments:
#             if len(fragment_pair) == 2:
#                 # If already broken into two parts, add it to the result
#                 result_fragments.append(fragment_pair)
#             elif len(fragment_pair) == 1:
#                 # If only one fragment, submit it for parallel recursion
#                 futures.append(executor.submit(recursive_bond_break, fragment_pair[0]))

#         # Collect the results from the parallel processes
#         for future in as_completed(futures):
#             result_fragments.extend(future.result())
        
#     result_fragments_list = [item for sublist in result_fragments for item in sublist]

#     return result_fragments_list

# def calculate_mol_weight(mol):
#     if mol is None:
#         return np.nan
    
#     return ExactMolWt(mol)
# def get_smiles(mol):
#     if mol is None:
#         return np.nan
    
#     return Chem.MolToSmiles(mol)
# import multiprocessing as mp
# def get_smiles_parallel(mol_list, num_processes=8):
#     with mp.Pool(processes=num_processes) as pool:
#         results = pool.map(get_smiles, mol_list)
#     return results
# def get_mw_parallel(mol_list, num_processes=8):
#     with mp.Pool(processes=num_processes) as pool:
#         results = pool.map(calculate_mol_weight, mol_list)
#     return results





# def prepare_molecule(smiles):
#     """Sanitize molecule to avoid Kekulization errors."""
#     mol = Chem.MolFromSmiles(smiles)
#     Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE)
#     return mol