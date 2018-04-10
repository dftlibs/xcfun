#.rst:
#
# Enables code coverage and defines the list of appropriate flags for the GNU
# compiler.
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--f90bindings=<F90BINDINGS> Test Fortran 90 bindings [default: ON]."
#     - "--pybindings Enable Python interface [default: OFF]."
#   define:
#     - "'-DTEST_Fortran_BINDINGS={0}'.format(arguments['--f90bindings'])"
#     - "'-DENABLE_PYTHON_INTERFACE={0}'.format(arguments['--pybindings'])"

option_with_print(
  NAME
    TEST_Fortran_BINDINGS
  MESSAGE
    "Test Fortran 90 bindings"
  DEFAULT
    ON
  )

option_with_print(
  NAME
    ENABLE_PYTHON_INTERFACE
  MESSAGE
    "Enable Python interface"
  DEFAULT
    OFF
  )

if(TEST_Fortran_BINDINGS)
  enable_language(Fortran)
  include(FortranCInterface)
  FortranCInterface_VERIFY(CXX)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
  add_subdirectory(fortran)
endif()

if(ENABLE_PYTHON_INTERFACE)
  add_subdirectory(python)
endif()
