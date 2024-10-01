from setuptools import setup, find_packages

VERSION = '0.0.2'
DESCRIPTION = 'Spectral denoising and denoising search'
LONG_DESCRIPTION = 'A package that allows to denoise MS/MS spectra based on molecular formula information and intensity modeling. Denoise search integrates spectral desnoing into identity search process'
REQUIREMENTS = [i.strip() for i in open("requirements.txt").readlines()]
# Setting up
setup(
    name="",
    version=VERSION,
    author="Fanzhou Kong",
    author_email="<fzkong@ucdavis.edu>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=REQUIREMENTS,
    keywords=['python', 'metabolomics', 'MSMS denoise'],
    entry_points = {
        "console_scripts":[
            'spectral_denoiisng_status = spectral_denoising:verify_status',
        ],
    },

)
