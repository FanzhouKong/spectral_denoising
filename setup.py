from setuptools import setup, find_packages

VERSION = '0.2.7'
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
    long_description= long_description,
    long_description_content_type='text/markdown',

)
