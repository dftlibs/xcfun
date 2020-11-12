# -*- coding: utf-8 -*-
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import os
import subprocess
import sys
import pathlib
import shutil
import re

# -- Project information -----------------------------------------------------

project = 'XCFun'
copyright = '2020, Ulf Ekström and contributors'
author = 'Ulf Ekström and contributors'

# The full version, including alpha/beta/rc tags
release = '2.1.1'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.mathjax', 'breathe', 'recommonmark']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The master toctree document.
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    html_theme = 'default'
else:
    html_theme = 'sphinx_rtd_theme'

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = '_static/logo2.png'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

breathe_projects = {'XCFun': '_build/xml'}
breathe_default_project = 'XCFun'


class Xlator(dict):
    """ All-in-one multiple-string-substitution class """
    def _make_regex(self):
        """ Build re object based on the keys of the current dictionary """
        return re.compile("|".join(map(re.escape, self.keys())))

    def __call__(self, match):
        """ Handler invoked for each regex match """
        return self[match.group(0)]

    def xlat(self, text):
        """ Translate text, returns the modified text. """
        return self._make_regex().sub(self, text)


def run_doxygen(folder):
    """Run the doxygen make command in the designated folder"""

    pathlib.Path("_build").mkdir(exist_ok=True)
    subs = Xlator({
        "@PROJECT_SOURCE_DIR@": str(pathlib.Path(__file__).resolve().parents[1]),
        "@PERL_EXECUTABLE@": str(shutil.which("perl")),
        "@DOXYGEN_DOT_PATH@": str(shutil.which("dot")),
    })
    doxy_in = (pathlib.Path(__file__).parent / "Doxyfile.in").resolve()
    doxy_out = pathlib.Path("_build/Doxyfile").resolve()
    with doxy_in.open("r") as f, doxy_out.open("w") as g:
        doxy = f.read()
        g.write(subs.xlat(doxy))

    try:
        retcode = subprocess.call(f"cd {folder}; {shutil.which('doxygen')}", shell=True)
        if retcode < 0:
            sys.stderr.write(f"doxygen terminated by signal {-retcode}")
    except OSError as e:
        sys.stderr.write(f"doxygen execution failed: {e}")


def setup(app):
    # symlink CHANGELOG.md
    confpy = pathlib.Path(__file__)
    changelog = confpy.parent / "changelog.md"
    if not changelog.exists():
        changelog.symlink_to(confpy.parents[1] / "CHANGELOG.md")
    run_doxygen('_build')
