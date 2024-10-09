import pandas as pd 
import re
import os
from . import spectral_operations as so
import re
import shutil
import numpy as np
from tqdm import tqdm
from fuzzywuzzy import fuzz
import json
def read_msp(file_path):
    
    """
    Reads the MSP files into the pandas dataframe, and sort/remove zero intensity ions in MS/MS spectra.

    Args:
        file_path (str): target path path for the MSP file.
    Returns:
        pd.DataFrame: DataFrame containing the MS/MS spectra information
    """
    
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
        # Save the last spectrum
        if spectrum:
            spectra.append(spectrum)
    df = pd.DataFrame(spectra)
    df['peaks'] = [so.sort_spectrum(so.remove_zero_ions(np.array(peak))) for peak in df['peaks']]
    for column in df.columns:
        if column != 'peaks':  # Skip 'peaks' column
            df[column] = pd.to_numeric(df[column], errors='ignore')
        # elif column == 'peaks':
        #     peak = 
    return df
def write_to_msp(df, file_path, msms_col = 'peaks'):
    
    """
    Pair function of read_msp.
    Exports a pandas DataFrame to an MSP file.

    Args:
        df (pd.DataFrame): DataFrame containing spectrum information. Should have columns for 'name', 'peaks', and other metadata.
        file_path (str): Destination path for the MSP file.
    Returns:
        None
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
def save_df(df, save_path):
    """
    Pair function of save_df.

    Save a DataFrame contaning MS/MS spectra to a CSV file, converting any columns containing 2D numpy arrays to string format.

    Args:
        df (pandas.DataFrame): The DataFrame to be saved.
        save_path (str): The file path where the DataFrame should be saved. If the path does not end with '.csv', it will be appended automatically.
    Returns:
        None
    Notes:
        - This function identifies columns in the DataFrame that contain 2D numpy arrays with a second dimension of size 2.
        - These identified columns are converted to string format before saving to the CSV file.
        - The function uses tqdm to display a progress bar while processing the rows of the DataFrame.
    """

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
    """
    Pair function of write_df.
    Reads a CSV file into a DataFrame, processes specific columns based on a pattern check, 
    and MS/MS in string format to 2-D numpy array (string is used to avoid storage issue in csv files).
    
    Args:
        path (str): The file path to the CSV file.
    Returns:
        pandas.DataFrame: The processed DataFrame with specific columns converted.
    Raises:
        FileNotFoundError: If the file at the specified path does not exist.
        pd.errors.EmptyDataError: If the CSV file is empty.
        pd.errors.ParserError: If the CSV file contains parsing errors.
    Notes:
        - The function assumes that the first row of the CSV file contains the column headers.
        - The `check_pattern` function is used to determine which columns to process.
        - The `so.str_to_arr` function is used to convert the values in the selected columns.
    """
    
    df = pd.read_csv(path)
    print('done read in df')
    cols = []
    for c in df.columns:
        if check_pattern(df.iloc[0][c]):
            cols.append(c)
    for col in cols:
        df[col] = [so.str_to_arr(y[col]) for x,y in df.iterrows()]
    return(df)

def sanitize_df_cols(df):

    ''' 
    Sanitize the column names of a DataFrame by removing the 'reference' prefix and replacing 'msms' with 'peaks'. Mainly used for compatibility with the data files from previous versions.
    
    Args:
        df (pandas.DataFrame): the dataframe to sanitize
    Returns:
        df (pandas.DataFrame): the dataframe with columns are sanitized for spectral denoising/denoising search
    '''
    
    df.columns = df.columns.str.replace('reference_', '', regex=True)
    df.columns = df.columns.str.replace('msms', 'peaks')
    return df
def check_pattern(input_string):
    """
    Helper function for read_df.
    Regular expression to match pairs of floats in standard or scientific notation separated by a tab, with each pair on a new line

    Args:
        input_string (str): input string to check for the pattern
    Returns:
        bool: True if the pattern is found, False otherwise
    """
    
    if isinstance(input_string, str):
        if '\t' in input_string:
            return True
    return False