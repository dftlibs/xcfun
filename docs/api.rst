.. _api:

*****************************************
XCFun's application programming interface
*****************************************

The library is written in C++, but can also be directly used in a C or
Fortran project through its application programming interface.
The C interface is exposed described in the ``api/xcfun.h``, while the
Fortran interface is described in the module file ``api/xcfun.f90``.
This documentation describes the C API. The Fortran API is written as a wrapper
to the C API and has the same behavior.

Types and type definitions
++++++++++++++++++++++++++

.. doxygenstruct:: XCFunctional

.. doxygentypedef:: xcfun_t

Functions
+++++++++

.. doxygenfunction:: xcfun_version

.. doxygenfunction:: xcfun_splash

.. doxygenfunction:: xcfun_authors

.. doxygenfunction:: xcfun_test

.. doxygenfunction:: xcfun_is_compatible_library

.. doxygenfunction:: xcfun_which_vars

.. doxygenfunction:: xcfun_which_mode

.. doxygenfunction:: xcfun_enumerate_parameters

.. doxygenfunction:: xcfun_enumerate_aliases

.. doxygenfunction:: xcfun_describe_short

.. doxygenfunction:: xcfun_describe_long

.. doxygenfunction:: xcfun_new

.. doxygenfunction:: xcfun_delete

.. doxygenfunction:: xcfun_set

.. doxygenfunction:: xcfun_get

.. doxygenfunction:: xcfun_is_gga

.. doxygenfunction:: xcfun_is_metagga

.. doxygenfunction:: xcfun_eval_setup

.. doxygenfunction:: xcfun_user_eval_setup

.. doxygenfunction:: xcfun_input_length

.. doxygenfunction:: xcfun_output_length

.. doxygenfunction:: xcfun_eval

.. doxygenfunction:: xcfun_eval_vec

Enumerations
++++++++++++

.. doxygenenum:: xcfun_mode

.. doxygenenum:: xcfun_vars

Preprocessor definitions and global variables
+++++++++++++++++++++++++++++++++++++++++++++

.. doxygendefine:: XCFUN_API_VERSION 

.. doxygendefine:: XCFUN_MAX_ORDER 

.. doxygenvariable:: XCFUN_TINY_DENSITY 

.. doxygenvariable:: XC_EORDER 

.. doxygenvariable:: XC_EVARS

.. doxygenvariable:: XC_EMODE
