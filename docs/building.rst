.. _building:

Building XCFun
==============

.. _dependencies:

Dependencies
------------

- A C++ compiler compliant with the C++11 standard. `See here
  <https://en.cppreference.com/w/cpp/compiler_support#cpp11>`_ for a list of
  compatible compilers.
- The `CMake <https://cmake.org>`_ build system generator. Version 3.11 or later
  is required. To install a recent version of CMake locally::

    $ CMAKE_VERSION=3.14.7
    $ target_path=$HOME/Deps/cmake/$CMAKE_VERSION
    $ cmake_url="https://cmake.org/files/v${CMAKE_VERSION%.*}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz"
    $ mkdir -p "$target_path"
    $ curl -Ls "$cmake_url" | tar -xz -C "$target_path" --strip-components=1
    $ export PATH=$HOME/Deps/cmake/$CMAKE_VERSION/bin${PATH:+:$PATH}

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

To compile the standalone examples:

- A Fortran compiler with complete ``iso_c_binding`` support.
- A C compiler compliant with the C99 standard.

To compile the Python bindings:

- Python 3.6+ and its development libraries and headers.
- `pybind11 <https://pybind11.readthedocs.io>`_. This will be automatically
  downloaded if not available.
- Other dependencies listed in the `requirements.txt` or `environment.yml` files.

To compile the documentation:

- `Doxygen <http://doxygen.nl/>`_
- `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`_
- The `Breathe <https://breathe.readthedocs.io>`_ Sphinx extension.
- The ``recommonmark`` Sphinx extension.

.. _confbuildtest:

Configuring, building, testing
------------------------------

1. Clone the repository from GitHub or download a tarball with the sources.
2. Configure::

     $ cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=<install-prefix>

   We also provide a Python script as front-end to CMake, see :ref:`compile-options`.

3. Build::

     $ cd build
     $ make

4. Test::

     $ ctest

5. Install::

     $ make install

**Congratulations**, you are all set to use XCFun! Read on for details on :ref:`using`.

.. _compile-options:

Compilation options
-------------------

A Python script called ``setup`` is made available as a front-end to CMake. The basic configuration command::

  $ cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=<install-prefix>

translates to the following invocation of the ``setup`` script::

  $ python setup --prefix=<install-prefix>

The script's options mirror exactly the options you can set by directly using CMake.

- ``--cxx`` / ``CMAKE_CXX_COMPILER``. The C++ compiler to use to compile the library.
- ``--type`` / ``CMAKE_BUILD_TYPE``. Any of the build types recognized by
  CMake, *i.e.* debug, release, and so forth.
- ``<build-dir>`` / ``-B<build-dir>``. The location of the build folder.
- ``--xcmaxorder`` / ``XCFUN_MAX_ORDER``. Maximum derivative order, defaults to 6.
- ``--pybindings`` / ``XCFUN_PYTHON_INTERFACE``. Enable compilation of Python
  bindings, defaults to ``OFF``.
- ``--static`` / ``BUILD_SHARED_LIBS``. Compile only the static library,
  defaults to ``OFF``, building the shared library only.
- ``ENABLE_TESTALL``. Whether to compile unit tests. ``ON`` by default. To
  toggle it ``OFF`` when using the ``setup`` script use
  ``--cmake-options="-DENABLE_TESTALL=OFF"``.

.. _building-docs:

Building the documentation
--------------------------

To build the documentation::

  $ cd docs
  $ make html

or::

  $ sphinx-build docs _build -t html

Bumping versions
----------------

To bump a version you should edit the ``cmake/custom/xcfun.cmake``,
``src/version_info.hpp``, and ``docs/conf.py`` files.
