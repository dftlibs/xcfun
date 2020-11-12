# Change Log

## [Version 2.1.1] - 2020-11-12

### Changed

- Linux and macOS continuous integration testing is run on GitHub actions. See [PR #145](https://github.com/dftlibs/xcfun/pull/145)

### Fixed

- We polished the installation of header files, CMake target export files, and Python module. These are especially relevant for Conda packaging XCFun. See [PR #143](https://github.com/dftlibs/xcfun/pull/143)
- A numerical issue with SCAN functionals and small density gradients was fixed by James Furness (@JFurness1). See [issue #144](https://github.com/dftlibs/xcfun/issues/144) reported by Xing Zhang (@fishjojo) and subsequent [PR #146](https://github.com/dftlibs/xcfun/pull/146) for the fix.

## [Version 2.1.0] - 2020-09-18

- Many new functionals in the SCAN family have been added. Thanks to James
  Furness for the contribution.
  See [PR #140](https://github.com/dftlibs/xcfun/pull/140)
- The library is now available both as a Spack and a Conda package.
- The library can now be _natively_ compiled on Linux, macOS, and Windows.

### Changed

- **BREAKING** CMake >= 3.14 is required to configure the code.

## [Version 2.0.2] - 2020-07-15

### Fixed

- VWN3 functional has been fixed for the spin-polarized case. It previously gave wrong results
  when alpha and beta densities differed. Thanks to Zhenyu Zhu for reporting the problem
  and also suggesting the solution.
  See [PR #134](https://github.com/dftlibs/xcfun/pull/134) and 
  [issue #132](https://github.com/dftlibs/xcfun/issues/132).

## [Version 2.0.1] - 2020-05-06

### Fixed

- We removed the `DEBUG_POSTFIX` property from the properties on the `xcfun`
  target. This was leading to build failures when using the library through
  CMake `FetchContent` with mixed release/debug mode.

## [Version 2.0.0] - 2020-04-14

### Changed

- **BREAKING** The build system will only produce a shared (default) or static
  library. Compilation of the static library can be requested by setting
  `BUILD_SHARED_LIBS` to `OFF`.
- macOS CI testing was moved to Azure Pipelines.
- The dependency on pybind11 was bumped to v2.5.0

### Fixed

- We corrected a number of wrinkles in the handling of symbol visibility in the
  shared library.

## [Version 2.0.0a7] - 2020-04-10

### Fixed

- Address warnings from compilers. Fix #90.

## [Version 2.0.0a6] - 2020-02-23

### Fixed

- Compilation with GCC 5.4.0.

## [Version 2.0.0a5] - 2020-02-20

### Fixed

- Handling of 64-bit integers in the Fortran interface.

## [Version 2.0.0a4] - 2020-02-02

### Fixed

- The API function `xcfun_get` accepts a single in-out `double` parameter. It
  was erroneously declared to accept an array of `double`-s instead.

## [Version 2.0.0a3] - 2020-01-31

We have introduced a number of breaking changes, motivated by the need to
modernize the library. See the [migration guide](https://xcfun.readthedocs.io/en/latest/migration.html).

### Added

- Up-to-date API documentation generated with [Doxygen], [breathe], and [Sphinx].
- Up-to-date documentation on how to build and develop XCFun.
- Up-to-date documentation on how to use XCFun in your code.
- API functions `xcfun_which_vars` and `xcfun_which_mode`.
- A full example, based on CMake as build system generator, showing how to use
  the library from a C++ host. Thanks @stigrj!
- A full example, based on CMake as build system generator, showing how to use
  the library from a C host.
- A full example, based on CMake as build system generator, showing how to use
  the library from a Fortran host.

### Changed

- **BREAKING** All API functions are uniformly namespaced with the `xcfun_` prefix.
- **BREAKING** The Fortran interface has been completely rewritten using
  `iso_c_binding`: the library can now be compiled without the use of neither a
  C nor a Fortran compiler. :confetti_ball:
- **BREAKING** CMake option `XCFun_XC_MAX_ORDER` renamed to `XCFUN_MAX_ORDER`. New default value of 6.
- **BREAKING** CMake option `XCFun_ENABLE_PYTHON_INTERFACE` renamed to `XCFUN_PYTHON_INTERFACE`.

### Deprecated

### Removed

- **BREAKING** API functions `xc_serialize`, `xc_deserialize`, `xc_set_fromstring`, and `xc_derivative_index`.
- **BREAKING** The CMake options `ENABLE_FC_SUPPORT` and `ENABLE_64BIT_INTEGERS`.

### Fixed

### Security

## [Version 2.0.0a2] - 2020-01-21

## [Version 2.0.0a1] - 2019-12-15

### Added

- A user-friendly API function to set up functional evaluation `xc_user_eval_setup`. Thanks @ilfreddy.

### Changed

- **BREAKING** A compiler compliant with the C++11 (or later) standard is required.
- **BREAKING** CMake >= 3.11 is required to configure the code. 
- **BREAKING** The Python bindings are now generated using [pybind11] instead of
  SWIG. The dependency will be fetched at configuration time if not found on
  your system.
- **BREAKING** The Fortran interface is no longer build with the code, but
  shipped as a separate file to be compiled within your own Fortran code.

[Unreleased]: https://github.com/dftlibs/xcfun/compare/v2.1.1...HEAD
[Version 2.1.1]: https://github.com/dftlibs/xcfun/compare/v2.1.0...v2.1.1
[Version 2.1.0]: https://github.com/dftlibs/xcfun/compare/v2.0.2...v2.1.0
[Version 2.0.2]: https://github.com/dftlibs/xcfun/compare/v2.0.1...v2.0.2
[Version 2.0.1]: https://github.com/dftlibs/xcfun/compare/v2.0.0...v2.0.1
[Version 2.0.0]: https://github.com/dftlibs/xcfun/compare/v2.0.0a7...v2.0.0
[Version 2.0.0a7]: https://github.com/dftlibs/xcfun/compare/v2.0.0a6...v2.0.0a7
[Version 2.0.0a6]: https://github.com/dftlibs/xcfun/compare/v2.0.0a5...v2.0.0a6
[Version 2.0.0a5]: https://github.com/dftlibs/xcfun/compare/v2.0.0a4...v2.0.0a5
[Version 2.0.0a4]: https://github.com/dftlibs/xcfun/compare/v2.0.0a3...v2.0.0a4
[Version 2.0.0a3]: https://github.com/dftlibs/xcfun/compare/v2.0.0a2...v2.0.0a3
[Version 2.0.0a2]: https://github.com/dftlibs/xcfun/compare/v2.0.0a1...v2.0.0a2
[Version 2.0.0a1]: https://github.com/dftlibs/xcfun/releases/tag/v2.0.0a1

[GitHub]: https://github.com/dftlibs/xcfun
[pybind11]: https://pybind11.readthedocs.io
[Doxygen]: http://doxygen.nl/
[breathe]: https://breathe.readthedocs.io/en/latest/
[Sphinx]: https://www.sphinx-doc.org/en/master/
