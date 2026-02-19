"""Sphinx configuration for Molecular Modeling Toolkit documentation."""

import os
import sys

# Add source directory to path so autodoc can find modules
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../source'))

project = 'Molecular Modeling Toolkit'
copyright = '2025, Molecular Modeling Development Team'
author = 'Molecular Modeling Development Team'
version = '0.1.0'
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'ase': ('https://wiki.fysik.dtu.dk/ase/', None),
}

# Default language
language = 'ja'

# Theme
html_theme = 'alabaster'
html_theme_options = {
    'description': 'ASE-based molecular modeling with Voronoi optimization',
    'github_button': False,
}

html_static_path = ['_static']
templates_path = ['_templates']

# Master doc (each language has its own index)
master_doc = 'ja/index'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Autodoc settings
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}
