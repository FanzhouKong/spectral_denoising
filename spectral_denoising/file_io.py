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
    if os.path.exists(file_path)== False:
        raise FileNotFoundError(f"File not found: {file_path}")
        return ()
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
            try:
                df[column] = pd.to_numeric(df[column], errors='raise')
            except:
                pass
    df = standardize_col(df) 
    return df
def write_to_msp(df, file_path, msms_col = 'peaks', normalize = False):
    
    """
    Pair function of read_msp.
    Exports a pandas DataFrame to an MSP file.

    Args:
        df (pd.DataFrame): DataFrame containing spectrum information. Should have columns for 'name', 'peaks', and other metadata.
        file_path (str): Destination path for the MSP file.
    Returns:
        None
    """
    if normalize == True:
        df[msms_col] = [so.normalize_spectrum(peak) for peak in df[msms_col]]
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
def read_df(path, keep_ms1_only = False):
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
    
    print('done read in df...')
    for col in df.columns:
        if check_pattern(df[col].iloc[0]):
            df[col] = [so.str_to_arr(y[col]) for x,y in df.iterrows()]
    df =  standardize_col(df)



    if keep_ms1_only == False:
        df.dropna(subset=['peaks'], inplace=True)
    if ':' in df.iloc[0]['peaks']:
        df['peaks']=[so.msdial_to_array(row['peaks']) for index, row in df.iterrows() if row['peaks'] == row['peaks']]
    df.reset_index(drop=True, inplace=True)
    return(df)

from .constant import standard_mapping
def standardize_col(df):
    """
    Standardizes column names in the given DataFrame based on a provided mapping. Help to read in and processing files with MS Dial generated msp files.
    
    Args:
    df (pd.DataFrame): The DataFrame whose column names need to be standardized.
    
    standard_mapping (dict): A dictionary where keys are common variations of the name, 
                              and values are the standard name.
    
    Returns:
    pd.DataFrame: DataFrame with standardized column names.
    """
    
    # Create a mapping for case-insensitive column names
    new_columns = []
    for col in df.columns:
        # Convert the column name to lowercase
        col_lower = col.lower()
        if col_lower != 'reference_precursor_mz':
            col_lower = col_lower.replace('reference_', '')
        # Map the column name to the standard one if found in the standard mapping
        standardized_col = standard_mapping.get(col_lower)
        if standardized_col is not None:
            new_columns.append(standardized_col)
        else:
            new_columns.append(col_lower)
    # Assign the new standardized columns back to the DataFrame
    df.columns = new_columns
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
def export_denoising_searches(results, save_dir, top_n = 10):
    """
    Pair function of import_denoising_searches.
    Exports the results of a denoising search to a JSON file.

    Args:
        results (list): The list of results from a denoising search.
        save_path (str): The file path where the results should be saved. If the path does not end with '.json', it will be appended automatically.
    Returns:
        None
    """
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    for i in range(len(results)):
        if results[i].empty or len(results[i])==0:
            continue
        else:
            temp = results[i].head(top_n)
            pmz_temp = temp.iloc[0]['precursor_mz']
            write_to_msp(temp, os.path.join(save_dir, f"denoising_search_{i}_{pmz_temp:0.4f}.msp"), msms_col='query_peaks_denoised')
