# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
import os
# sys.path.insert(0, os.path.abspath("../../"))
# sys.path.insert(0, os.path.abspath("../../spectral_denoising"))
project = 'Spectral denoising'
copyright = '2024, Fanzhou Kong'
author = 'Fanzhou Kong'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", 
              "sphinx.ext.doctest", 
              "sphinx.ext.autosummary", 
              "sphinx.ext.viewcode", 
              "sphinx.ext.githubpages",
              "sphinxcontrib.email",
              "numpydoc"]

numpydoc_show_class_members = False
autoclass_content = 'both'

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]