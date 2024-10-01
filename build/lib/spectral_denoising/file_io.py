import pandas as pd 
import re
import os
import spectral_denoising.spectral_operations as so
import re
import shutil
import numpy as np
from tqdm import tqdm
from fuzzywuzzy import fuzz
import json
def read_msp(file_path):
    spectra = []
    spectrum = {}

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            
            # Handle metadata
            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip().lower()
                value = value.strip()
                
                if key == 'name':
                    # Save current spectrum and start a new one
                    if spectrum:
                        spectra.append(spectrum)
                    spectrum = {'name': value, 'peaks': []}
                else:
                    spectrum[key] = value
            
            # Handle peak data (assumed to start with a number)
            elif line[0].isdigit():
                
                peaks = line.split()
                m_z = float(peaks[0])
                intensity = float(peaks[1])
                spectrum['peaks'].append((([m_z, intensity])))
        # spectrum['peaks'] = np.array(spectrum['peaks'])
        # spectrum['peaks'] = so.sort_spectra(so.remove_zero_ions(spectrum['peaks']))
        # Save the last spectrum
        if spectrum:
            spectra.append(spectrum)
    df = pd.DataFrame(spectra)
    df['peaks'] = [so.sort_spectra(so.remove_zero_ions(np.array(peak))) for peak in df['peaks']]
    for column in df.columns:
        if column != 'peaks':  # Skip 'peaks' column
            df[column] = pd.to_numeric(df[column], errors='ignore')
        # elif column == 'peaks':
        #     peak = 
    return df
def write_to_msp(df, file_path, msms_col = 'peaks'):
    """
    Exports a pandas DataFrame to an MSP file.

    Parameters:
    df (pd.DataFrame): DataFrame containing spectrum information.
                       Should have columns for 'name', 'peaks', and other metadata.
    file_path (str): Destination path for the MSP file.
    """
    with open(file_path, 'w') as f:
        for _, row in df.iterrows():
            # Write the name of the spectrum
            if isinstance(row[msms_col], float):
                continue
            if 'name' in df.columns:
                f.write(f"Name: {row['name']}\n")
            
            # Write other metadata if available
            for col in df.columns:
                if col not in ['name', msms_col] and 'peak' not in col:
                    f.write(f"{col.capitalize()}: {row[col]}\n")
            
            # Write the peaks (assuming each peak is a tuple of (m/z, intensity))
            f.write(f"Num Peaks: {len(row[msms_col])}\n")
            for mz, intensity in row[msms_col]:
                f.write(f"{mz} {intensity}\n")
            
            # Separate spectra by an empty line
            f.write("\n")
def check_pattern(input_string):
    # Regular expression to match pairs of floats in standard or scientific notation separated by a tab,
    # with each pair on a new line
    if isinstance(input_string, str):
        if '\t' in input_string:
            return True
    return False
def save_df(df, save_path):
    data = df.copy()
    cols = []
    for c in df.columns:
        if isinstance(df.iloc[0][c], np.ndarray):
            if np.shape(df.iloc[0][c])[1]==2:
                cols.append(c)
    print(cols)
    if save_path.endswith('.csv') == False:
        save_path = save_path+'.csv'
    for col in cols:
        specs = []
        for index, row in tqdm(data.iterrows(), total = len(data)):
            specs.append(so.arr_to_str(row[col]))
        data[col]=specs
    data.to_csv(save_path, index = False)
def read_df(path):
    df = pd.read_csv(path)
    print('done read in df')
    cols = []
    for c in df.columns:
        if check_pattern(df.iloc[0][c]):
            cols.append(c)
    for col in cols:
        df[col] = [so.str_to_arr(y[col]) for x,y in df.iterrows()]
    return(df)
import pandas as pd
def sanitize_df_cols(df):
    # Remove 'reference_' prefix from all column names
    df.columns = df.columns.str.replace('^reference_', '', regex=True)
    df.columns = df.columns.str.replace('msms', 'peaks')
    return df