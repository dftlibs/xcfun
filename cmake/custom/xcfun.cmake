#.rst:
#
# Sets general XCFun options
#
# Variables modified::
#
#   XCFun_XC_MAX_ORDER -- Maximum order of derivatives of the exchange-correlation kernel
#   ENABLE_PYTHON_INTERFACE -- Whether to enable the Python interface
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--xcmaxorder=<XCFun_XC_MAX_ORDER> An integer greater than 3 [default: 3]."
#     - "--pybindings Enable Python interface [default: OFF]."
#   define:
#     - "'-DXCFun_XC_MAX_ORDER=\"{0}\"'.format(arguments['--xcmaxorder'])"
#     - "'-DENABLE_PYTHON_INTERFACE={0}'.format(arguments['--pybindings'])"

option_with_default(XCFun_XC_MAX_ORDER "Maximum order of derivatives of the exchange-correlation kernel" 3)
# Make sure user selected a valuer larger than 2
if(DEFINED XCFun_XC_MAX_ORDER AND XCFun_XC_MAX_ORDER LESS 3)
  message(STATUS "${XCFun_XC_MAX_ORDER} not a valid value for maximum order of XC kernel derivatives! Resetting to its default value 3")
  set(XCFun_XC_MAX_ORDER 3 CACHE STRING "Maximum order of derivatives of the exchange-correlation kernel" FORCE)
endif()

set(PROJECT_VERSION 2.0.0)
set(PROJECT_VERSION_MAJOR 2)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)

# Hardcode to share, rather than use CMAKE_INSTALL_DATAROOTDIR as the latter
# might resolve to a place not recognized by CMake
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")

if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
endif()
file(TO_NATIVE_PATH "${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}/xcfun" PYMOD_INSTALL_FULLDIR)

option_with_print(ENABLE_PYTHON_INTERFACE "Enable Python interface" OFF)

if(ENABLE_FC_SUPPORT)
  enable_language(Fortran)
  include(FortranCInterface)
  FortranCInterface_VERIFY(CXX)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
  include(${CMAKE_CURRENT_LIST_DIR}/compilers/FortranFlags.cmake)
  include(${CMAKE_CURRENT_LIST_DIR}/int64.cmake)
endif()

if(ENABLE_PYTHON_INTERFACE)
  add_subdirectory(python)
endif()

add_subdirectory(api)
