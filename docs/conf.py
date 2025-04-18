# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# #Add modules for automatic documentation

# %%

# make shure that the package directory is in the sys path
import os, sys
from pathlib import Path


# Note that on the github workflow the build is executed in the doc folder,
# Thus if that is the case, we need to go up the foldertree

curfolder = os.path.abspath(".")
if "docs" in curfolder:
    # when executing in docs folder
    basefolder = Path(curfolder).parents[0]

else:
    # when executing in basefolder
    basefolder = curfolder


sys.path.insert(0, str(basefolder))
# sys.path.insert(0, os.path.join(str(basefolder), "metobs_toolkit"))


# Test importing on github workflow
if "/runner/" in os.getcwd():
    print("ASSUME SERVER BUILD OF DOCUMENTATION")
    print(f"sys.path: {sys.path}")
    import metobs_toolkit

    print(
        f" the version of the toolkit to create docs for: {metobs_toolkit.__version__}"
    )
    # sys.path.insert(0, '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit')


# The toolkit must be imported when testing and building the documentation
# locally. However this is overkill for RTD service, so only import it for
# local builds
if "/home/thoverga" in str(basefolder):
    import metobs_toolkit

logofile = os.path.join(basefolder, "docs", "logo_wide_1280x640.jpeg")


# %%


# -- Project information -----------------------------------------------------

project = "metobs_toolkit"
copyright = "2023, Thomas Vergauwen"
author = "Thomas Vergauwen"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",  # Autodocument functions
    "sphinx_rtd_theme",  # Use the read the docs theme
    "sphinx.ext.viewcode",  # Button to go to source code
    "sphinx_copybutton",  # Copy button (for examples etc)
    "sphinx.ext.napoleon",  # To convert Numpydocstring to readable format
    "sphinx.ext.autosummary",  # Create neat summary tables
    "myst_parser",  # for including md files (readme)
    "sphinx.ext.autosectionlabel",  # for cross linking
    "nbsphinx",  # to render the notebook examples in the doc
    "matplotlib.sphinxext.plot_directive",  # embedding figures in the docstrings examples
]

# -- General configuration ------------------------------------------------------------

# Add any paths that contain templates here, relative to this directory.


templates_path = ["_templates"]
autosummary_generate = True  # Turn on sphinx.ext.autosummary

# Specify how to render the following file formats:
source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}
source_encoding = "utf-8"
master_doc = "index"  # The master toctree document.


# When building the doc, sphinx will try to import all the depending packages,
# This is needed because of the plot examples that are rendered in the docs!
# THus do not mock any package
# autodoc_mock_imports = [
#     "ee",
#     "pytz",
#     "matplotlib",
#     "numpy",
#     "geopandas",
#     "pandas",
#     "pyproj",
#     "shapely",
#     "cartopy",
#     "branca",
#     "geemap",
#     "folium",
#     "mpl_toolkits",
#     "scipy",
# ]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "paper/paper.md",
]

add_function_parentheses = False
add_module_names = False
show_authors = False  # section and module author directives will not be shown
todo_include_todos = False  # Do not show TODOs in docs


# Make sure the target is unique
autosectionlabel_prefix_document = True
# -- Options for HTML output -------------------------------------------------


html_theme = "pydata_sphinx_theme"
html_title = "MetObs Toolkit documentation"
html_short_title = "MetObs Toolkit documentation"
html_logo = "logo_small.svg"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sphinx = True
html_show_copyright = True
htmlhelp_basename = "MetObs toolkit"  # Output file base name for HTML help builder.
html_use_smartypants = True
html_show_sourcelink = True


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/vergauwenthomas/MetObs_toolkit",
            "icon": "fab fa-github-square fa-xl",
        },
    ],
    # "navbar_center": ["version-switcher", "navbar-nav"],
    # "switcher": {
    #     # The json url must be a full path operationally !!!
    #     "json_url": "https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/docs/_static/custom.css",  # this file contains a dict of all versions to show
    #     "version_match": f"v{version}",  # currently being browsed
    # },
}


# Try to remove the white/dark theme switch (because it is bugging)
# html_context = {
#    "default_mode": "light" #light or dark theme
# }


# html_theme_options["navbar_end"] = ["navbar-icon-links"]

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = f"MetObs Toolkit {metobs_toolkit.__version__} documentation"

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = "MetObs Toolkit documentation"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = "logo_small.svg"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = "_static/logo/favicon.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
# html_css_files = [
#     "custom.css",
# ]

# =============================================================================
# nbsphinx settings
# =============================================================================

# isue is that some example notebooks require the GEE authentication,
# which fails when documentation is build online !!

# since the execution of the notebooks is part of the development pipeline,
# we can skip (in general) this step and assume the notebooks have output.

# but, since this package is under active development, it is handy that the
# notbooks are executed only when building locally !!

if ("/runner/" in os.getcwd()) | ("readthedocs.org" in os.getcwd()):
    print("ASSUME SERVER BUILD OF DOCUMENTATION")

    # The notebooks are executed BEFORE the documentation is build!
    nbsphinx_execute = "never"  # never, always or auto
else:
    print("ASSUME LOCAL BUILD OF DOCUMENTATION")
    nbsphinx_allow_errors = True  # for developping
    nbsphinx_execute = "never"  # never, always or auto

# =============================================================================
# Matplotlib include settings
# =============================================================================

# See: https://matplotlib.org/stable/api/sphinxext_plot_directive_api.html

plot_include_source = True
plot_html_show_source_link = False
plot_formats = ["png"]  # no need for highres and pdf versions

plot_pre_code = "import numpy as np \nfrom matplotlib import pyplot as plt\nplt.rcParams['figure.autolayout'] = True"
