import random
import numpy as np
from . import spectral_operations as so
def generate_noise(pmz, lamda, n = 100):
    """
    Generate synthetic electronic noise for spectral data.
    
    Parameters:
        pmz (float): The upper bound for the mass range.

        lamda (float): The lambda parameter for the Poisson distribution, which serves as both mean and standard deviation of the distribution.
        
        n (int, optional): The number of random noise ions to generate. Defaults to 100.
    Returns:
        np.array: A synthetic spectrum with electronic noise.
    """

    if int(n)!= n:
        n = np.int64(np.ceil(n))
    else:
        n = n
    # Generate a random variable from a uniform distribution in the range [a, b]
    mass = [random.uniform(50, pmz) for _ in range(n)]

    # size specifies the number of random variates to generate.

    # Generating Poisson-distributed random variables
    intensity = np.random.poisson(lam=lamda, size=n)
    intensity = intensity/100
    return(so.pack_spectrum(mass, intensity))
def add_noise(msms, noise):
    """
    Add noise to a mass spectrum and process the resulting spectrum.
    This function takes a mass spectrum and a noise spectrum, standardizes the mass spectrum,
    adds the noise to it, normalizes the resulting spectrum, and sorts it.

    Args:
        msms (np.ndarray): The mass spectrum to which noise will be added.

        noise (np.ndarray): The noise spectrum to be added to the mass spectrum.

    Returns:
        np.ndarray: The processed mass spectrum after adding noise, normalization, and sorting.
    
    Notes:
        - The noise spectrum is generated with intensity as ralatie measure (from 0-1)
        - Thus, the mass spectrum is standardized using the standardize_spectrum function.

    """

    msms = so.standardize_spectrum(msms)
    msms_c = so.add_spectra(msms, noise)
    return(so.sort_spectrum(so.normalize_spectrum(msms_c)) )
def generate_chemical_noise(pmz, lamda, polarity,formula_db,n = 100):
    
    """
    Generate chemical noise for a given mass-to-charge ratio (m/z) and other parameters.
    The m/z of the chemical noise is taken from a database of all true possible mass values. 
    The detailes about this database can be found paper: LibGen: Generating High Quality Spectral Libraries of Natural Products for EAD-, UVPD-, and HCD-High Resolution Mass Spectrometers

    Args:
        pmz (float): The target mass-to-charge ratio (m/z) value.
        
        lamda (float): The lambda parameter for the Poisson distribution used to generate intensities, which serves as both mean and standard deviation of the distribution.

        polarity (str): The polarity of the adduct, either '+' or '-'.

        formula_db (pandas.DataFrame): A DataFrame containing a column 'mass' with possible mass values.
        
        n (int, optional): The number of noise peaks to generate. Default is 100.

    Returns:
        np.array: A synthetic spectrum with chemical noise.

    Raises:
        ValueError: If the polarity is not '+' or '-'.
    """

    mass_e =  -0.00054858026
    if polarity =='+':
        coe = 1
    elif polarity =='-':
        coe = -1
    else:
        print('cannot determine adduct polarity!')
        return()
    if int(n)!= n:
        n = np.int64(np.ceil(n))
    else:
        n = n
    all_possible_mass = np.array(formula_db['mass'])
    idx_left, idx_right = all_possible_mass.searchsorted([50,pmz ])
    all_allowed_mass = all_possible_mass[idx_left:idx_right]
    if idx_right-idx_left <n:
        n = idx_right-idx_left
    # Generate a random variable from a uniform distribution in the range [a, b]
    mass = np.random.choice(all_allowed_mass, size=n, replace=False)
    mass = mass+coe*mass_e
    # mass = [random.uniform(50, pmz) for _ in range(n)]

    # size specifies the number of random variates to generate.

    # Generating Poisson-distributed random variables
    intensity = np.random.poisson(lam=lamda, size=n)
    intensity = intensity/100
    return(so.pack_spectrum(mass, intensity))
