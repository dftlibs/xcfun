.. _api:

*****************************************
XCFun's application programming interface
*****************************************

The library is written in C++, but can also be directly used in a C or
Fortran project through its application programming interface.
The C interface is exposed described in the ``api/xcfun.h``, while the
Fortran interface is described in the module file ``api/xcfun.F90``.
This documentation describes the C API. The Fortran API is written as a wrapper
to the C API and has the same behavior.

Functions
+++++++++

.. doxygenfunction:: xcfun_version
   :project: XCFun

.. doxygenfunction:: xcfun_splash
   :project: XCFun

.. doxygenfunction:: xc_new_functional_not_macro
   :project: XCFun

.. doxygenfunction:: xc_free_functional
   :project: XCFun

.. doxygenfunction:: xc_serialize
   :project: XCFun

.. doxygenfunction:: xc_enumerate_parameters
   :project: XCFun

.. doxygenfunction:: xc_enumerate_aliases
   :project: XCFun

.. doxygenfunction:: xc_set
   :project: XCFun

.. doxygenfunction:: xc_get
   :project: XCFun

.. doxygenfunction:: xc_describe_short
   :project: XCFun

.. doxygenfunction:: xc_describe_long
   :project: XCFun

.. doxygenfunction:: xc_is_gga
   :project: XCFun

.. doxygenfunction:: xc_is_metagga
   :project: XCFun

.. doxygenfunction:: xc_set_fromstring
   :project: XCFun

.. doxygenfunction:: xc_user_eval_setup
   :project: XCFun

.. doxygenfunction:: xc_eval_setup
   :project: XCFun

.. doxygenfunction:: xc_eval
   :project: XCFun

.. doxygenfunction:: xc_eval_vec
   :project: XCFun

.. doxygenfunction:: xc_derivative_index
   :project: XCFun

Enumerations
++++++++++++

.. doxygenenum:: xc_mode
   :project: XCFun

.. doxygenenum:: xc_vars
   :project: XCFun

Preprocessor definitions
++++++++++++++++++++++++

.. doxygendefine:: XC_MAX_ORDER 
   :project: XCFun

.. doxygendefine:: XC_TINY_DENSITY 
   :project: XCFun

.. doxygendefine:: XC_NO_REGULARIZATION 
   :project: XCFun

.. doxygendefine:: XC_EORDER 
   :project: XCFun

.. doxygendefine:: XC_EVARS
   :project: XCFun

.. doxygendefine:: XC_EMODE
   :project: XCFun
