[![Build Status](https://travis-ci.org/robertodr/xcfun.svg?branch=master)](https://travis-ci.org/robertodr/xcfun)
[![Build status](https://ci.appveyor.com/api/projects/status/78xtl84v5gufjs9l?svg=true)](https://ci.appveyor.com/project/miroi/xcfun)

# Arbitrary order Exchange-Correlation functional library

Copyright Ulf Ekstrom <uekstrom@gmail.com> and contributors 2009-2010.
See http://dftlibs.org/xcfun/ for more information.

The main interface is in `include/xcfun.h`
(or `fortran/xc_fun_module.f90` for Fortran bindings).

## Copying

The library is licensed under the LGPL license version 3, see
COPYING.LESSER for more information.


## Configuration

Check that `XC_MAX_ORDER` is defined to the highest order derivatives
you need (and not higher) in `src/config.h~.  Using a too large value
for `XC_MAX_ORDER` makes compilation slow and the generated code huge.


## Building a debug/development version

Edit the Makefile that matches your compiler
(Makefile.gcc or Makefile.intel or ...)
to set CXX (C++ compiler) and flags and run

    $ make -f Makefile.gcc

(or using the corresponding Makefile)
This will create the library file `lib/libxcfun.so`


## Building an optimized version

Edit the Makefile that matches your compiler
(Makefile.gcc or Makefile.intel or ...)
and add `-DNDEBUG` to the compiler flags. Add optimization
compiler options. Make sure your compiler performs inlining
(-O3 with gcc).

# Primary test environments

## Continuous integration builds

The Travis CI builds are triggered only when pushing to the `master` branch.

- Ubuntu 12.04 LTS 64-bit with Python 2.7.3 and CMake 3.3.2
  this is the environment offered by [Travis CI](https://travis-ci.org) pulling
  in various PPA. The following compilers are used, both in release and debug:

  1. GCC 4.6
  2. GCC 4.7
  3. GCC 4.8
  4. GCC 4.9
  5. GCC 5.1
  6. Clang 3.5 and GFortran 4.6
  7. Clang 3.6 and GFortran 4.6
  8. Clang 3.7 and GFortran 4.6
  9. Clang 3.8 and GFortran 4.6

- Mac OS X 10.9.5 with Python 2.7.10 and CMake 3.2.3
  this is the environment offered by [Travis CI](https://travis-ci.org)
  The following compilers are used, both in release and debug:

  1. XCode 6.4 with Clang and GFortran 5.2
  2. XCode 6.4 with GCC 5.2
  3. XCode 7.0 with Clang and GFortran 5.2
  4. XCode 7.0 with GCC 5.2
