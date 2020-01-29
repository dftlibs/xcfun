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

Functions
+++++++++

.. doxygenfunction:: xcfun_version
   :project: XCFun

.. doxygenfunction:: xcfun_splash
   :project: XCFun

.. doxygenfunction:: xcfun_authors
   :project: XCFun

.. doxygenfunction:: xcfun_test
   :project: XCFun

.. doxygenfunction:: xcfun_is_compatible_library
   :project: XCFun

.. doxygenfunction:: xcfun_which_vars
   :project: XCFun

.. doxygenfunction:: xcfun_which_mode
   :project: XCFun

.. doxygenfunction:: xcfun_enumerate_parameters
   :project: XCFun

.. doxygenfunction:: xcfun_enumerate_aliases
   :project: XCFun

.. doxygenfunction:: xcfun_describe_short
   :project: XCFun

.. doxygenfunction:: xcfun_describe_long
   :project: XCFun

.. doxygenfunction:: xcfun_new
   :project: XCFun

.. doxygenfunction:: xcfun_delete
   :project: XCFun

.. doxygenfunction:: xcfun_set
   :project: XCFun

.. doxygenfunction:: xcfun_get
   :project: XCFun

.. doxygenfunction:: xcfun_is_gga
   :project: XCFun

.. doxygenfunction:: xcfun_is_metagga
   :project: XCFun

.. doxygenfunction:: xcfun_eval_setup
   :project: XCFun

.. doxygenfunction:: xcfun_user_eval_setup
   :project: XCFun

.. doxygenfunction:: xcfun_input_length
   :project: XCFun

.. doxygenfunction:: xcfun_output_length
   :project: XCFun

.. doxygenfunction:: xcfun_eval
   :project: XCFun

.. doxygenfunction:: xcfun_eval_vec
   :project: XCFun

Enumerations
++++++++++++

.. doxygenenum:: xcfun_mode
   :project: XCFun

.. doxygenenum:: xcfun_vars
   :project: XCFun

Preprocessor definitions and global variables
+++++++++++++++++++++++++++++++++++++++++++++

.. doxygendefine:: XCFUN_API_VERSION 
   :project: XCFun

.. doxygendefine:: XCFUN_MAX_ORDER 
   :project: XCFun

.. doxygenvariable:: XCFUN_TINY_DENSITY 
   :project: XCFun

.. doxygenvariable:: XC_EORDER 
   :project: XCFun

.. doxygenvariable:: XC_EVARS
   :project: XCFun

.. doxygenvariable:: XC_EMODE
   :project: XCFun
