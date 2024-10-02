import random
import numpy as np
from .spectral_operations import standardize_spectrum, sort_spectrum, add_spectra,normalize_spectrum, pack_spectrum
def generate_noise(pmz, lamda, n = 100):
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
    return(pack_spectrum(mass, intensity))
def add_noise(msms, noise):
    msms = standardize_spectrum(msms)
    msms_c = add_spectra(msms, noise)
    return(sort_spectrum(normalize_spectrum(msms_c)) )
def generate_chemical_noise(pmz, lamda, polarity,formula_db,n = 100):
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
    return(pack_spectrum(mass, intensity))
