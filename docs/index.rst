.. XCFun documentation master file, created by
   sphinx-quickstart on Fri Nov 18 14:51:49 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to XCFun documentation
==============================

XCFun is a library of exchange-correlation (XC) functionals to be used in
electronic-structure theory codes.
XCFun follows a unique implementation strategy which enables the computation of
derivatives of the XC functional kernel up to arbitrary order. It does so by
relying on forward-mode automatic differentiation.

Given a new XC functional kernel, its implementation with all its derivatives
only requires to write code for the undifferentiated kernel.
This implementation strategy is very powerful and allows:

1. Faster implementation of new functionals: you write the kernel, the compiler
   does the rest.
2. Introduction of new variables, for example current densities, in the
   parametrization of new or existing XC kernels.
3. Testing the numerical stability of XC kernels.

Contents:

.. toctree::
   :maxdepth: 2

   xcfun_install.rst
   migration.rst
   api.rst
   functionals.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

