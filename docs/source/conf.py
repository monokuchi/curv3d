# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Geometry'
copyright = '2023'
author = 'Kenneth Jao'

release = '0.1'
version = '0.1.0'

# -- General configuration

html_static_path = ["_static"]
html_js_files = ["js/mathjax-config.js"]

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.graphviz',
    'sphinxcontrib.bibtex',
    'breathe',
]

breathe_projects = {}
breathe_projects[project] = "../build/xml"
breathe_default_project = project

bibtex_bibfiles = ['refs.bib']

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,
}

# -- Options for EPUB output
epub_show_urls = 'footnote'
