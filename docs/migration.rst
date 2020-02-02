.. _migration:

*******************************************************
Migrating to the new application programmers' interface
*******************************************************

This is a short guide to migrating to the new application programmers' interface
(API) and build system for XCFun. We assume that you have successfully built and
tested XCFun and installed it to ``<install-prefix>``.
The layout of the install tree will be as follows:

.. code-block:: bash

   <install-prefix>/
                   ├── include
                   │   └── XCFun
                   │       ├── config.hpp
                   │       ├── densvars.hpp
                   │       ├── functional.hpp
                   │       ├── functionals
                   │       ├── specmath.hpp
                   │       ├── XCFunctional.hpp
                   │       ├── XCFunExport.h
                   │       ├── xcfun.f90
                   │       ├── xcfun.h
                   │       └── xcint.hpp
                   ├── lib
                   │   └── python
                   │       └── xcfun
                   ├── lib64
                   │   ├── libxcfun.a
                   │   ├── libxcfun.so -> libxcfun.so.2
                   │   └── libxcfun.so.2
                   └── share
                       └── cmake
                           └── XCFun

C/C++ host programs
-------------------

Types and function signatures are in the header file ``xcfun.h``.

In your source code, apply the following changes:

- Remove any of the calls to the functions that have been removed from the API.
  **Open an issue** if these functions are essential to your workflow and you would
  like them to be reinstated.
- Replace ``xc_functional`` with ``xcfun_t *``.
- Replace ``xc_new_functional`` with ``xcfun_new``.
- Replace ``xc_enumerate_parameters`` with ``xcfun_enumerate_parameters``.
- Replace ``xc_enumerate_aliases`` with ``xcfun_enumerate_aliases``.
- Replace ``xc_set`` with ``xcfun_set``.
- Replace ``xc_get`` with ``xcfun_get``.
- Replace ``xc_describe_short`` with ``xcfun_describe_short``.
- Replace ``xc_describe_long`` with ``xcfun_describe_long``.
- Replace ``xc_is_gga`` with ``xcfun_is_gga``.
- Replace ``xc_is_metagga`` with ``xcfun_is_metagga``.
- Replace ``xc_eval_setup`` with ``xcfun_eval_setup``.
- Replace ``xc_user_eval_setup`` with ``xcfun_user_eval_setup``.
- Replace ``xc_input_length`` with ``xcfun_input_length``.
- Replace ``xc_output_length`` with ``xcfun_output_length``.
- Replace ``xc_eval`` with ``xcfun_eval``.
- Replace ``xc_eval_vec`` with ``xcfun_eval_vec``.

Fortran host programs
---------------------

The Fortran/C interoperability layer for types and function signatures is in the source file ``xcfun.f90``.

In your source code, apply the following changes:

- Use the intrinsic ``iso_c_binding`` module: ``use, intrinsic :: iso_c_binding``.
- Remove any of the calls to the functions that have been removed from the API.
  **Open an issue** if these functions are essential to your workflow and you would
  like them to be reinstated.
- You should call the intrinsic ``trim`` on functions returning strings: ``xcfun_version``,
  ``xcfun_splash``, ``xcfun_authors``, ``xcfun_enumerate_paramters``,
  ``xcfun_enumerate_aliases``, ``xcfun_describe_short``,
  ``xcfun_describe_long``.
- Replace the type for the ``xc_functional`` object (now ``xcfun_t *``) from ``integer`` to ``type(c_ptr)``.
- Replace ``xc_new_functional`` with ``xcfun_new``.
- Replace ``xc_enumerate_parameters`` with ``xcfun_enumerate_parameters``.
- Replace ``xc_enumerate_aliases`` with ``xcfun_enumerate_aliases``.
- Replace ``xc_set`` with ``xcfun_set``.
- Replace ``xc_get`` with ``xcfun_get``.
- Replace ``xc_describe_short`` with ``xcfun_describe_short``.
- Replace ``xc_describe_long`` with ``xcfun_describe_long``.
- Replace ``xc_is_gga`` with ``xcfun_is_gga``.
- Replace ``xc_is_metagga`` with ``xcfun_is_metagga``.
- Replace ``xc_eval_setup`` with ``xcfun_eval_setup``.
- Replace ``xc_user_eval_setup`` with ``xcfun_user_eval_setup``.
- Replace ``xc_input_length`` with ``xcfun_input_length``.
- Replace ``xc_output_length`` with ``xcfun_output_length``.
- Replace ``xc_eval`` with ``xcfun_eval``.
- Replace ``xc_eval_vec`` with ``xcfun_eval_vec``.
