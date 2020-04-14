#.rst:
#
# Sets general XCFun options
#
# Variables modified::
#
#   XCFUN_MAX_ORDER -- Maximum order of derivatives of the exchange-correlation kernel
#   XCFUN_PYTHON_INTERFACE -- Whether to enable the Python interface
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--xcmaxorder=<XCFUN_MAX_ORDER> An integer greater than 3 [default: 6]."
#     - "--pybindings Enable Python interface [default: OFF]."
#   define:
#     - "'-DXCFUN_MAX_ORDER=\"{0}\"'.format(arguments['--xcmaxorder'])"
#     - "'-DXCFUN_PYTHON_INTERFACE={0}'.format(arguments['--pybindings'])"

option_with_default(XCFUN_MAX_ORDER "Maximum order of derivatives of the exchange-correlation kernel" 6)
# Make sure user selected a valuer larger than 3
if(DEFINED XCFUN_MAX_ORDER AND XCFUN_MAX_ORDER LESS 3)
  message(STATUS "${XCFUN_MAX_ORDER} not a valid value for maximum order of XC kernel derivatives! Resetting to its default value 6")
  set(XCFUN_MAX_ORDER 6 CACHE STRING "Maximum order of derivatives of the exchange-correlation kernel" FORCE)
endif()

set(PROJECT_VERSION 2.0.0)
set(PROJECT_VERSION_MAJOR 2)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKEDIR CMake)
else()
  set(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
endif()
set(CMAKECONFIG_INSTALL_DIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")

add_subdirectory(${PROJECT_SOURCE_DIR}/api)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

option_with_print(XCFUN_PYTHON_INTERFACE "Enable Python interface" OFF)

if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
endif()
file(TO_NATIVE_PATH "lib/${PYMOD_INSTALL_LIBDIR}/xcfun" PYMOD_INSTALL_FULLDIR)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})

if(XCFUN_PYTHON_INTERFACE)
  include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_pybind11.cmake)
  add_subdirectory(python)
endif()
