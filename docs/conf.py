"""Sphinx configuration for Molecular Modeling Toolkit documentation."""

import os
import platform
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

# Language — override with SPHINX_LANG env var (e.g. SPHINX_LANG=en)
language = os.environ.get('SPHINX_LANG', 'ja')

# Theme
html_theme = 'alabaster'
html_theme_options = {
    'description': 'ASE-based molecular modeling with Voronoi optimization',
    'github_button': False,
}

html_static_path = ['_static']
templates_path = ['_templates']

# Master doc — follows language setting
master_doc = f'{language}/index'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# --- LaTeX / PDF output settings ---
latex_engine = 'lualatex'
if language == 'ja':
    latex_documents = [
        ('ja/index', 'MolecularModelingToolkit_ja.tex',
         'Molecular Modeling Toolkit ドキュメント',
         author, 'manual'),
    ]
    latex_docclass = {
        'manual': 'ltjsbook',
        'howto': 'ltjsreport',
    }
else:
    latex_documents = [
        ('en/index', 'MolecularModelingToolkit_en.tex',
         'Molecular Modeling Toolkit Documentation',
         author, 'manual'),
    ]
    latex_docclass = {}
latex_use_xindy = False
latex_domain_indices = False

# Platform-aware font configuration for LuaLaTeX
_is_macos = platform.system() == 'Darwin'

if _is_macos:
    _fontpkg = r'''
\usepackage{luatexja-fontspec}
\setmainfont[BoldFont={Times New Roman Bold}]{Times New Roman}
\setsansfont{Helvetica}
\setmonofont{Menlo}
\setmainjfont[ItalicFont={Hiragino Mincho ProN}]{Hiragino Mincho ProN}
\setsansjfont[ItalicFont={Hiragino Sans}]{Hiragino Sans}
% Tell luatexja these symbols are NOT CJK — use the western font instead
\ltjdefcharrange{9}{"21A9, "21AA, "21BB, "2192, "2190,
  "2716, "2726, "25B6, "25C0, "25B2, "25BC}
\ltjsetparameter{jacharrange={-9}}
% Set western fonts that have these symbols
\newfontfamily\symbolfont{Apple Symbols}
\newfontfamily\menlosymbol{Menlo}
\usepackage{newunicodechar}
\newunicodechar{↩}{{\symbolfont ↩}}
\newunicodechar{↪}{{\symbolfont ↪}}
\newunicodechar{↻}{{\symbolfont ↻}}
\newunicodechar{✖}{{\symbolfont ✖}}
\newunicodechar{▶}{{\symbolfont ▶}}
\newunicodechar{◀}{{\symbolfont ◀}}
\newunicodechar{▲}{{\symbolfont ▲}}
\newunicodechar{▼}{{\symbolfont ▼}}
\newunicodechar{→}{{\symbolfont →}}
\newunicodechar{←}{{\symbolfont ←}}
\newunicodechar{✦}{{\menlosymbol ✦}}
'''
    _preamble = r'''
\usepackage{pict2e}
\usepackage{ellipse}
% Suppress small-caps italic warning (Times New Roman has no scit shape)
\DeclareFontShape{TU}{TimesNewRoman(0)}{m}{scit}{<->ssub*TimesNewRoman(0)/m/it}{}
\DeclareFontShape{TU}{TimesNewRoman(0)}{m}{sc}{<->ssub*TimesNewRoman(0)/m/n}{}
'''
else:
    # Linux (CI): use Noto CJK fonts (include Latin glyphs)
    _fontpkg = r'''
\usepackage{luatexja-fontspec}
\setmainfont{Noto Serif CJK JP}
\setsansfont{Noto Sans CJK JP}
\setmonofont{Noto Sans Mono CJK JP}
\setmainjfont{Noto Serif CJK JP}
\setsansjfont{Noto Sans CJK JP}
\ltjdefcharrange{9}{"21A9, "21AA, "21BB, "2192, "2190,
  "2716, "2726, "25B6, "25C0, "25B2, "25BC}
\ltjsetparameter{jacharrange={-9}}
\usepackage{newunicodechar}
\newunicodechar{↩}{{\textrm{↩}}}
\newunicodechar{↪}{{\textrm{↪}}}
\newunicodechar{↻}{{\textrm{↻}}}
\newunicodechar{✖}{{\textrm{✖}}}
\newunicodechar{▶}{{\textrm{▶}}}
\newunicodechar{◀}{{\textrm{◀}}}
\newunicodechar{▲}{{\textrm{▲}}}
\newunicodechar{▼}{{\textrm{▼}}}
\newunicodechar{→}{{\textrm{→}}}
\newunicodechar{←}{{\textrm{←}}}
\newunicodechar{✦}{{\textrm{✦}}}
'''
    _preamble = r'''
\usepackage{pict2e}
'''

latex_elements = {
    'polyglossia': '',
    'fontpkg': _fontpkg,
    'preamble': _preamble,
    'geometry': r'\usepackage[margin=25mm,footskip=15mm]{geometry}',
    'printindex': '',
}

# Autodoc settings
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}
