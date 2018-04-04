[![Build Status](https://travis-ci.org/dftlibs/xcfun.svg?branch=master)](https://travis-ci.org/dftlibs/xcfun)
[![GitHub license](https://img.shields.io/github/license/dftlibs/xcfun.svg?style=flat-square)](https://github.com/dftlibs/xcfun/blob/master/LICENSE.md)

[![GitHub issues](https://img.shields.io/github/issues/dftlibs/xcfun.svg?style=flat-square)](https://github.com/dftlibs/xcfun/issues)
[![GitHub forks](https://img.shields.io/github/forks/dftlibs/xcfun.svg?style=flat-square)](https://github.com/dftlibs/xcfun/network)
[![GitHub stars](https://img.shields.io/github/stars/dftlibs/xcfun.svg?style=flat-square)](https://github.com/dftlibs/xcfun/stargazers)

# Arbitrary order Exchange-Correlation functional library

Copyright [Ulf Ekström] and [contributors] 2009-2018. See
http://dftlibs.org/xcfun/ for more information.
The documentation is available at https://xcfun.readthedocs.io

The main interface is in `include/xcfun.h` (or `fortran/xc_fun_module.f90` for
Fortran bindings).

## License

XCFun is licensed under version 2.0 of the Mozilla Public License ([MPLv2.0]),
see [LICENSE.md].

## Contribution guide

- All patches are contributed through fork-pull request mechanism.
- We do not require any formal copyright assignment or contributor license
  agreement. Any contributions intentionally sent upstream are presumed to be
  offered under terms of the [MPLv2.0].
- Maintainers do not push directly to the repository.
- Maintainers do not review their own patches.
- Approval of one maintainer is sufficient for trivial patches.
  Trivial patches are typos and trivial bugfixes.
- For any patch that is not trivial, two maintainers need to sign off the patch.

## Configuration

Check that `XC_MAX_ORDER` is defined to the highest order derivatives you need
(and not higher) in `src/config.h`.
**Warning**: Using a too large value for `XC_MAX_ORDER` makes compilation slow
and the generated code huge.

## Building a debug/development version

Edit the Makefile that matches your compiler (Makefile.gcc or Makefile.intel or
...) to set CXX (C++ compiler) and flags and run

    $ make -f Makefile.gcc

(or using the corresponding Makefile) This will create the library file
`lib/libxcfun.so`

## Building an optimized version

Edit the Makefile that matches your compiler (Makefile.gcc or Makefile.intel or
...) and add `-DNDEBUG` to the compiler flags.
Add optimization compiler options. Make sure your compiler performs inlining
(-O3 with gcc).

[Ulf Ekström]: mailto:uekstrom@gmail.com
[contributors]: https://github.com/dftlibs/xcfun/blob/master/AUTHORS.md
[MPLv2.0]: https://www.mozilla.org/en-US/MPL/2.0/
[LICENSE.md]: https://github.com/dftlibs/xcfun/blob/master/LICENSE.md
