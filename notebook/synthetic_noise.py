import pandas as pd
from toolsets.file_io import read_df, save_df
import os
from tqdm import tqdm
import numpy as np
import toolsets.denoising_related_functions as drf
import toolsets.spectra_operations as so
from toolsets.search import string_search
master_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/dilution_series_data'
formula_db = pd.read_csv('/Users/fanzhoukong/Documents/GitHub/Libgen_stable/db/formulaDB_sorted.csv')
data = read_df(os.path.join(master_dir, 'dilution_mapped_unprocessed.csv'))
from itertools import product
chemical_lamda = 50
electric_lamda = 5
chemcial_abundant = [0.1,0.2,0.5]
electric_abundant = [2,10,100]
combinations = product(chemcial_abundant, electric_abundant)

# Convert to list and print each combination (optional)
for combination in tqdm(combinations):
    #first always chemical
    msms_conts_all = []
    entropy_raw_all = []
    entropy_denoised_all = []
    msms_denoised_all = []
    chemical_ratio = combination[0]
    electri_ratio = combination[1]
    data_current = data.copy()

    for index, row in data_current.iterrows()
        noise_chemical = so.generate_chemical_noise(row['reference_precursor_mz'],chemical_lamda, row['reference_adduct'][-1], formula_db=formula_db, n = np.int64(np.ceil(chemical_ratio*num_peaks)))
        noise_electric = so.generate_noise(row['reference_precursor_mz'], electric_lamda, n = np.int64(np.ceil(electri_ratio*num_peaks)))

    for f in tqdm(mrm_validated['processed_formulas'].unique()):
        data_temp = string_search(mrm_validated, 'processed_formulas', f)
        element_dict, element_count, element_mass, all_possible_mass, all_possible_candidate_formula = drf.get_all_subsets(f)
        all_formulas = [drf.dict_to_formula(x, element_dict) for x in all_possible_candidate_formula]
        allowed = np.array([drf.check_ratio(x) for x in all_formulas])
        all_element_mass = all_possible_mass[allowed]
        all_possible_candidate_formula =all_possible_candidate_formula[allowed]

        for index, row in data_temp.iterrows():
            diff_row = []
            entropy_raw_row = []
            entropy_denoised_row = []
            msms_conts_row = []
            msms_denoised_row = []
            num_peaks = so.num_peaks(row['msms'])
            for i in range(0, 100):
                noise_chemical = so.generate_chemical_noise(row['reference_precursor_mz'],chemical_lamda, row['reference_adduct'][-1], formula_db=formula_db, n = np.int64(np.ceil(chemical_ratio*num_peaks)))
                noise_electric = so.generate_noise(row['reference_precursor_mz'], electric_lamda, n = np.int64(np.ceil(electri_ratio*num_peaks)))
                msms_cont = so.add_noise(row['msms'], noise_chemical)
                msms_cont = so.add_noise(msms_cont, noise_electric)
                msms_d = drf.denoising_with_subset(msms_cont, row['reference_precursor_mz'], element_dict, all_possible_mass, all_possible_candidate_formula)
                entropy_raw = so.entropy_identity(msms_cont, row['library_peaks'], pmz = row['reference_precursor_mz'])
                entropy_raw_row.append(entropy_raw)
                entropy_denoised = so.entropy_identity(msms_d, row['library_peaks'], pmz = row['reference_precursor_mz'])
                entropy_denoised_row.append(entropy_denoised)
                diff_row.append(entropy_denoised-entropy_raw)
                msms_conts_row.append(msms_cont)
                msms_denoised_row.append(msms_d)
            optimal_idx = np.argmax(diff_row)
            msms_conts_all.append(msms_conts_row[optimal_idx])
            msms_denoised_all.append(msms_denoised_row[optimal_idx])
            entropy_raw_all.append(entropy_raw_row[optimal_idx])
            entropy_denoised_all.append(entropy_denoised_row[optimal_idx])
        result_df = pd.concat([result_df, data_temp], ignore_index=True)

    result_df['msms_cont']=msms_conts_all
    result_df['msms_cont_denoised']=msms_denoised_all
    result_df['entropy_cont_raw'] = entropy_raw_all
    result_df['entropy_cont_denoised'] = entropy_denoised_all
    file_name = 'synthetic_noise_sn_'+str(chemical_ratio)+'_'+str(electri_ratio)+'.csv'
    save_df(result_df, os.path.join(master_dir, file_name) )
