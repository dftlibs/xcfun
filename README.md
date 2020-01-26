[![Build Status](https://travis-ci.org/dftlibs/xcfun.svg?branch=master)](https://travis-ci.org/dftlibs/xcfun)
[![GitHub license](https://img.shields.io/github/license/dftlibs/xcfun.svg?style=flat-square)](https://github.com/dftlibs/xcfun/blob/master/LICENSE.md)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3576419.svg)](https://doi.org/10.5281/zenodo.3576419)

[![GitHub issues](https://img.shields.io/github/issues/dftlibs/xcfun.svg?style=flat-square)](https://github.com/dftlibs/xcfun/issues)
[![GitHub forks](https://img.shields.io/github/forks/dftlibs/xcfun.svg?style=flat-square)](https://github.com/dftlibs/xcfun/network)
[![GitHub stars](https://img.shields.io/github/stars/dftlibs/xcfun.svg?style=flat-square)](https://github.com/dftlibs/xcfun/stargazers)

# Arbitrary order Exchange-Correlation functional library

Copyright [Ulf Ekström] and [contributors] 2009-2018. See
http://dftlibs.org/xcfun/ for more information.
The documentation is available at https://xcfun.readthedocs.io

The main interface is in `api/xcfun.h` (or `api/xcfun.f90` for Fortran bindings).

## License

XCFun is licensed under version 2.0 of the Mozilla Public License ([MPLv2.0]),
see [LICENSE.md].

## Configuration

Check that `XCFUN_MAX_ORDER` is defined to the highest order derivatives you need
(and not higher) in `src/config.h`.
**Warning**: Using a too large value for `XCFUN_MAX_ORDER** makes compilation slow
and the generated code huge.

## Building the library

**TO BE WRITTEN***

[Ulf Ekström]: mailto:uekstrom@gmail.com
[contributors]: https://github.com/dftlibs/xcfun/blob/master/AUTHORS.md
[MPLv2.0]: https://www.mozilla.org/en-US/MPL/2.0/
[LICENSE.md]: https://github.com/dftlibs/xcfun/blob/master/LICENSE.md
