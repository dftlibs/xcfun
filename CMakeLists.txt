cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(XCFun
  LANGUAGES CXX C
  VERSION 2.0.0
  )

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# C99
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Layout according to the GNU standards
include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# Options handling utilities
include(CMakeParseArguments)
include(CMakeDependentOption)
# Macro for printing an option in a consistent manner
# Written by Lori A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard)
# Syntax: print_option(NAME <option to print> DEFAULT <was specified>)
macro(print_option)
  set(options)
  set(oneValueArgs NAME DEFAULT)
  set(multiValueArgs)
  cmake_parse_arguments(print_option
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
  if(NOT DEFINED ${print_option_NAME} OR "${${print_option_NAME}}" STREQUAL "")
    message(STATUS "Setting (unspecified) option ${print_option_NAME}: ${print_option_DEFAULT}")
  else()
    message(STATUS "Setting option ${print_option_NAME}: ${print_option_DEFAULT}")
  endif()
endmacro()
# Wraps an option with default ON/OFF. Adds nice messaging to option()
# Written by Lori A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard)
# Syntax: option_with_print(NAME <option name> MESSAGE <description> DEFAULT <default value>)
macro(option_with_print)
  set(options)
  set(oneValueArgs NAME MESSAGE DEFAULT)
  set(multiValueArgs)
  cmake_parse_arguments(option_with_print
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
  print_option(
    NAME ${option_with_print_NAME}
    DEFAULT ${option_with_print_DEFAULT}
    )
  option(${option_with_print_NAME} ${option_with_print_MESSAGE} ${option_with_print_DEFAULT})
endmacro()
# Wraps an option with a default other than ON/OFF and prints it
# Written by Lori A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard)
# NOTE: Can\'t combine with above b/c CMake handles ON/OFF options specially
# NOTE2: CMake variables are always defined so need to further check for if
#       they are the NULL string. This is also why we need the force
# Syntax: option_with_default(NAME <option name> MESSAGE <description> DEFAULT <default value>)
macro(option_with_default)
  set(options)
  set(oneValueArgs NAME MESSAGE DEFAULT)
  set(multiValueArgs)
  cmake_parse_arguments(option_with_default
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
  print_option(
    NAME ${option_with_default_NAME}
    DEFAULT "${option_with_default_DEFAULT}"
    )
  if(NOT DEFINED ${option_with_default_NAME} OR "${${option_with_default_NAME}}" STREQUAL "")
    set(${option_with_default_NAME} "${option_with_default_DEFAULT}" CACHE STRING ${option_with_default_MESSAGE} FORCE)
  endif()
endmacro()

if(NOT DEFINED CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type set to ${CMAKE_BUILD_TYPE}")

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
string(TOUPPER "${CMAKE_BUILD_TYPE}" cmake_build_type_toupper)

if(NOT cmake_build_type_tolower STREQUAL "debug" AND
   NOT cmake_build_type_tolower STREQUAL "release" AND
   NOT cmake_build_type_tolower STREQUAL "minsizerel" AND
   NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
    message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo, and MinSizeRel (case-insensitive).")
endif()

# List options
option_with_print(
  NAME
    ENABLE_64BIT_INTEGERS
  MESSAGE
    "Enable 64-bit integers"
  DEFAULT
    OFF
  )

option_with_default(
  NAME
    XCFun_XC_MAX_ORDER
  MESSAGE
    "Maximum order of derivatives of the exchange-correlation kernel"
  DEFAULT
    3
  )
# Make sure user selected a valuer larger than 2
if(DEFINED XCFun_XC_MAX_ORDER AND XCFun_XC_MAX_ORDER LESS 3)
  message(STATUS "${XCFun_XC_MAX_ORDER} not a valid value for maximum order of XC kernel derivatives! Resetting to its default value 3")
  set(XCFun_XC_MAX_ORDER 3 CACHE STRING "Maximum order of derivatives of the exchange-correlation kernel" FORCE)
endif()
# FIXME Localize use of this option
add_definitions("-DXC_MAX_ORDER=${XCFun_XC_MAX_ORDER}")

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/taylor
    ${CMAKE_CURRENT_SOURCE_DIR}/src/functionals
    )

set(BASE_SOURCES
    src/xcint.cpp
    src/xcfun.cpp
    src/fortran.c
    )

file(GLOB FUNCTIONAL_SOURCES src/functionals/*.cpp)

add_library(
    xcfun
    ${BASE_SOURCES}
    ${FUNCTIONAL_SOURCES}
    )
# Set .so version
if(BUILD_SHARED_LIBS)
 set_target_properties(xcfun PROPERTIES
   VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION} )
endif()

install(TARGETS xcfun DESTINATION lib${LIB_SUFFIX})
file(GLOB headers "include/*.h")
install(FILES ${headers} DESTINATION include)

# FIXME We'll take care of compilers and compiler flags later
#if(ENABLE_FORTRAN_INTERFACE)
#    include(FortranCompilers)
#endif()
#include(CCompilers)
#include(CXXCompilers)

option_with_print(
  NAME
    ENABLE_Fortran_INTERFACE
  MESSAGE
    "Enable Fortran interface"
  DEFAULT
    ON
  )
if(ENABLE_Fortran_INTERFACE)
  enable_language(Fortran)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
  add_library(
    xcfun_f90_bindings
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran/xcfun_module.F90
    )

  # Set .so version
  if(BUILD_SHARED_LIBS)
    set_target_properties(xcfun_f90_bindings PROPERTIES VERSION ${VERSION} SOVERSION ${SOVERSION} )
  endif()
  install(TARGETS xcfun_f90_bindings DESTINATION lib${LIB_SUFFIX})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/xcfun.mod DESTINATION include)

  if(ENABLE_64BIT_INTEGERS)
    add_definitions("-DXCFUN_INTEGER=integer*8")
    add_definitions("-DXCFUN_FORTRAN_INT=long long int")
  endif()

  add_executable(
    kernel_example
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran/kernel_example.F90
    )
  target_link_libraries(
    kernel_example
    xcfun_f90_bindings
    xcfun
    )
  set_property(TARGET kernel_example PROPERTY LINKER_LANGUAGE Fortran)
endif()

option_with_print(
  NAME
    ENABLE_PYTHON_INTERFACE
  MESSAGE
    "Enable Python interface"
  DEFAULT
    OFF
  )
if(ENABLE_PYTHON_INTERFACE)
  find_package(SWIG REQUIRED)
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  find_package(NumPy REQUIRED)

  set_property(TARGET xcfun PROPERTY POSITION_INDEPENDENT_CODE ON)

  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_SOURCE_DIR}/python")

  set(XCFUN_SWIG         "${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun_swig.i")
  set(XCFUN_SWIG_WRAP    "${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swigPYTHON_wrap.c")
  set(CMAKE_SWIG_OUTDIR   ${CMAKE_CURRENT_BINARY_DIR}/python)

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swigPYTHON_wrap.c
                     COMMAND "${SWIG_EXECUTABLE}"
                     ARGS "-python" "-I${CMAKE_CURRENT_SOURCE_DIR}/include" "-outdir" "${CMAKE_SWIG_OUTDIR}" "-o" ${XCFUN_SWIG_WRAP} ${XCFUN_SWIG}
                     DEPENDS xcfun ${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun_swig.i)

  set(XCFUN_PY       "${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun.py")
  set(XCFUN_SWIG_PY  "${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swig.py")
  set(XCFUN_PY_LIB   "${CMAKE_CURRENT_BINARY_DIR}/python/_xcfun_swig.so")

  set(SETUP_PY_IN  "${CMAKE_CURRENT_SOURCE_DIR}/python/setup.py.in")
  set(SETUP_PY     "${CMAKE_CURRENT_BINARY_DIR}/python/setup.py")

  configure_file(${SETUP_PY_IN} ${SETUP_PY})

  add_custom_command(OUTPUT ${XCFUN_PY_LIB}
                     COMMAND ${PYTHON_EXECUTABLE}
                     ARGS "setup.py" "build_ext" "--inplace"
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python
                     DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/python/setup.py ${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swigPYTHON_wrap.c)

  add_custom_target(target ALL DEPENDS ${XCFUN_PY_LIB})

  install(FILES "${XCFUN_PY}"           DESTINATION "python")
  install(FILES "${XCFUN_SWIG_PY}"      DESTINATION "python")
  install(FILES "${XCFUN_PY_LIB}"       DESTINATION "python")
endif()


option_with_print(
  NAME
    ENABLE_TESTALL
  MESSAGE
    "Enable compilation of testall"
  DEFAULT
    ON
  )
if(ENABLE_TESTALL)
  add_executable(
      testall
      test/testall.c
      )
  target_link_libraries(
      testall
      xcfun
      )
  add_test(testall ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testall)
endif()

add_test(kernel_example ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kernel_example)
if(ENABLE_PYTHON_INTERFACE)
  add_test(NAME python-interface  COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test/test.py")
endif()

include(CTest)
enable_testing()
