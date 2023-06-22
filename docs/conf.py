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
sys.path.insert(0, os.path.join(str(basefolder), "metobs_toolkit"))

try:
    import metobs_toolkit
except:
    print("NOT ABLE TO IMPORT THE TOOLKIT!!")
    pass

print(sys.path)
# %%


# -- Project information -----------------------------------------------------

project = "metobs_toolkit"
copyright = "2023, Thomas Vergauwen"
author = "Thomas Vergauwen"

# The full version, including alpha/beta/rc tags
release = "0.0.1"


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
    "sphinx.ext.autosectionlabel", #for cross linking
]


# Add any paths that contain templates here, relative to this directory.


templates_path = ["_templates"]
autosummary_generate = True  # Turn on sphinx.ext.autosummary

# When building the doc, sphinx will try to import all the depending packages,
# this is not needed and problematic when building the docs in a clean docker on gitlab.
# So specify which packages can be mocked

autodoc_mock_imports = [
    "ee",
    "pytz",
    "matplotlib",
    "numpy",
    "geopandas",
    "pandas",
    "pyproj",
    "shapely",
    "cartopy",
]


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

html_theme_options = {
    "analytics_id": "G-XXXXXXXXXX",  #  Provided by Google in your dashboard
    "analytics_anonymize_ip": False,
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    "vcs_pageview_mode": "",
    "style_nav_header_background": "white",
    # Toc options
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_static_path = ["_static"]
