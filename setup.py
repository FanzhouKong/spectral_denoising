from setuptools import setup, find_packages

VERSION = '0.0.8'
DESCRIPTION = 'Spectral denoising and denoising search'
REQUIREMENTS = [i.strip() for i in open("requirements.txt").readlines()]
with open("README.md", "r") as file:
    long_description = file.read()
# Setting up
setup(
    name="",
    version=VERSION,
    author="Fanzhou Kong",
    author_email="<fzkong@ucdavis.edu>",
    description=DESCRIPTION,
    packages=find_packages(),
    install_requires=REQUIREMENTS,
    keywords=['python', 'metabolomics', 'MSMS denoise'],
    entry_points = {
        "console_scripts":[
            'spectral_denoiisng_status = spectral_denoising:verify_status',
        ],
    },
    long_description= long_description,
    long_description_content_type='text/markdown',

)
