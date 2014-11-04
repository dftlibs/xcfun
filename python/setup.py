#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

import os.path
import numpy

from distutils.core import setup, Extension

xcfun_swig_module = Extension('_xcfun_swig',
                              sources=['xcfun_swig_wrap.c'],
                              include_dirs=[numpy.get_include()],
                              library_dirs=['../../build'],
                              libraries=['xcfun', 'stdc++']
                              )

setup (name = 'example',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [xcfun_swig_module],
       py_modules = ["xcfun_swig"],
       )
