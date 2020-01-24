.. _building:

Building XCFun
==============

Compile the library using one of the included makefiles (for example
``Makefile.gcc``). This will generate a library file ``lib/libxcfun.a`` which
can be statically linked to your application. Include files listing all
available functionals will be generated during the compilation of XCFun. These
files are ``include/xcfun_autogen.h`` (for the C interface), and
``fortran/xcfun_autogen.f90`` for the Fortran module. C or C++ programs that
uses XCFun should include the ``xcfun.h`` header file, while Fortran programs
should use the xcfun module defined in ``fortran/xcfun_module.f90``.

See also Compile Time Options.
