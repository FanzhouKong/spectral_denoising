============
Installation
============

Spectral denoising requires ``Python >= 3.8`` installed on your system. It has been tested on Windows, Linux, and macOS platforms.


Installing from PyPI
====================

To install the latest version of MS Entropy from PyPI, use the following command:

.. code-block:: bash

  pip install spectral-denoising

MS Entropy package also provides many useful functions like reading .mzML/.msp/.mgf/.lbm2 files, to keep the package lightweight, these functions are not installed by default. If you need to use these functions, you can install the package with the ``all`` extra:


Dependencies
============

If you install Spectral denoising from PyPI, the necessary dependencies will be installed automatically. However, if you are installing from the source, install these manually:

- ``chemparse==0.3.1``
- ``CIRpy==1.0.2``
- ``fuzzywuzzy==0.18.0``
- ``matplotlib==3.9.2``
- ``molmass==2021.6.18``
- ``ms_entropy==1.3.3``
- ``numexpr==2.10.1``
- ``numpy==2.1.1``
- ``pandas==2.2.3``
- ``plotly==5.24.1``
- ``PubChemPy==1.0.4``
- ``rdkit==2024.3.5``
- ``Requests==2.32.3``
- ``scikit_learn==1.5.2``
- ``scipy==1.14.1``
- ``seaborn==0.13.2``
- ``tqdm==4.66.5``


Installing from Source
======================

Spectral denoising is designed to work in python environment. Preferentially in IPython or Jupyter Notebook.

- Downlaod and compile

  .. code-block:: bash

    git clone https://github.com/FanzhouKong/spectral_denoising.git
    cd spectral_denoising
    pip install -r requirements.txt
    
  Then, as in Pure Python Mode, you can either copy the spectral_denoising folder to your project directory or include it in your PYTHONPATH environment variable.


