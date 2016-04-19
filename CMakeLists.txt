cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(XCFun)

# XCFun version
set(VERSION 2.0.0)
# Library version
set(SOVERSION 2.0.0)

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers" OFF)

# on some systems it is difficult or expensive to install Fortran
# and some projects that use XCFun, do not need the Fortran interface (example: XCint)
# so we make it possible to disable it
# by default we build the Fortran interface
option(ENABLE_FORTRAN_INTERFACE "Enable Fortran interface" ON)

# Python interface is disabled by default, can be switched on with DENABLE_PYTHON_INTERFACE=On
option(ENABLE_PYTHON_INTERFACE "Enable Python interface" OFF)

if(ENABLE_FORTRAN_INTERFACE)
    enable_language(CXX C Fortran)
else()
    enable_language(CXX C)
endif()

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

include(ConfigSafeGuards)

# set maximum derivative order
set(XC_MAX_ORDER "3" CACHE STRING "Maximum derivative order")
add_definitions(-DXC_MAX_ORDER=${XC_MAX_ORDER})

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
  VERSION ${VERSION} SOVERSION ${SOVERSION} )
endif(BUILD_SHARED_LIBS)

install(TARGETS xcfun DESTINATION lib${LIB_SUFFIX})
file(GLOB headers "include/*.h")
install(FILES ${headers} DESTINATION include)

add_executable(
    testall
    test/testall.c
    )
target_link_libraries(
    testall
    xcfun
    )

if(ENABLE_FORTRAN_INTERFACE)
    include(FortranCompilers)
endif()
include(CCompilers)
include(CXXCompilers)

if(ENABLE_FORTRAN_INTERFACE)
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

if (ENABLE_PYTHON_INTERFACE)
    FIND_PACKAGE(SWIG REQUIRED)
    FIND_PACKAGE(PythonInterp REQUIRED)
    FIND_PACKAGE(PythonLibs REQUIRED)
    FIND_PACKAGE(NumPy REQUIRED)

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

    INSTALL(FILES "${XCFUN_PY}"           DESTINATION "python")
    INSTALL(FILES "${XCFUN_SWIG_PY}"      DESTINATION "python")
    INSTALL(FILES "${XCFUN_PY_LIB}"       DESTINATION "python")

endif()

add_test(testall        ${CMAKE_CURRENT_BINARY_DIR}/testall)
add_test(kernel_example ${CMAKE_CURRENT_BINARY_DIR}/kernel_example)
if(ENABLE_PYTHON_INTERFACE)
    add_test(NAME python-interface  COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test/test.py")
endif()

include(CTest)
enable_testing()
