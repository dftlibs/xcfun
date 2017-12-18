[![Build Status](https://travis-ci.org/dftlibs/xcfun.svg?branch=master)](https://travis-ci.org/dftlibs/xcfun)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](http://www.gnu.org/licenses/lgpl-3.0)

# Arbitrary order Exchange-Correlation functional library

Copyright Ulf Ekström <uekstrom@gmail.com> and [contributors](https://github.com/dftlibs/xcfun/blob/master/AUTHORS.md) 2009-2017.
See http://dftlibs.org/xcfun/ for more information.

The main interface is in `include/xcfun.h`
(or `fortran/xc_fun_module.f90` for Fortran bindings).

## License

The project license is specified in [COPYING.md] and [COPYING.LESSER.md]

XCFun is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

In addition to the rights granted to you under the LGPL license
(described in [COPYING.md] and [COPYING.LESSER.md]) you may distribute the
compiled XCFun library statically linked to an application. It is
understood that doing so will prevent the end user from relinking the
application to a newer version of the XCFun library, but static
linking is explicitly allowed by this exception. This exception does
not free you from your obligation under the LGPL to provide the full
(possibly modified by you) source code of the XCFun library to anyone
who receives a copy of the library or an application linked to the
library.

XCFun is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

[COPYING.md]: https://github.com/dftlibs/xcfun/blob/master/COPYING.md 
[COPYING.LESSER.md]: https://github.com/dftlibs/xcfun/blob/master/COPYING.LESSER.md


## Documentation

- http://xcfun.readthedocs.io


## Contribution guide

- All patches are contributed through fork-pull request mechanism.
- We do not require any formal copyright assignment or contributor license
  agreement. Any contributions intentionally sent upstream are presumed to be
  offered under terms of the LGPL version 3.0.
- Maintainers do not push directly to the repository.
- Maintainers do not review their own patches.
- Approval of one maintainer is sufficient for trivial patches.
- Trivial patches are typos and trivial bugfixes.
- For any patch that is not trivial, two maintainers need to sign off the patch.
- **TODO**: We need to better define what "trivial" means.

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
