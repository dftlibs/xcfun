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
option(ENABLE_PYTHON_INTERFACE "Enable Python interface" ON)

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

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/taylor
    ${CMAKE_CURRENT_SOURCE_DIR}/src/functionals
    ${CMAKE_CURRENT_SOURCE_DIR}/python
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
        add_definitions(-DXCFUN_FORTRAN_INT='long long int')
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
    INCLUDE(${SWIG_USE_FILE})

    FIND_PACKAGE(PythonLibs)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/python)

    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_SOURCE_DIR}/python")

    SET(CMAKE_SWIG_FLAGS "-v")
    SET(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/python)

    SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun_swig.i PROPERTIES CPLUSPLUS OFF)

    SWIG_ADD_MODULE(xcfun_swig python ${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun_swig.i ${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swigPYTHON_wrap.c)
    SWIG_LINK_LIBRARIES(xcfun_swig ${PYTHON_LIBRARIES} xcfun)

    FIND_PROGRAM(PYTHON "python" )
# now for the part doing the python setup
    if (PYTHON)

        SET(PYTHON_DEST "${CMAKE_CURRENT_BINARY_DIR}/python")
        INSTALL(TARGETS _xcfun_swig DESTINATION ${PYTHON_DEST})
        INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun.py  DESTINATION ${PYTHON_DEST})

        set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/python/setup.py.in")
        set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/python/setup.py")
        set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/python")

        configure_file(${SETUP_PY_IN} ${SETUP_PY})

        add_custom_command(OUTPUT ${OUTPUT}
                           COMMAND ${PYTHON}
                           ARGS "setup.py build_ext --inplace"
                           DEPENDS ${DEPS})

        add_custom_target(target ALL DEPENDS ${OUTPUT})

        execute_process(COMMAND "${PYTHON} ${SETUP_PY}" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/python")

    endif()

endif()

include(Tests)
include(CTest)
enable_testing()
