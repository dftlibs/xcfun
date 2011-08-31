cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(XCFun CXX C Fortran)

# XCFun version
set(VERSION 2.0.0)
# Library version
set(SOVERSION 2.0.0)


set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

# set build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
                AND NOT cmake_build_type_tolower STREQUAL "release"
                AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
        message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
        Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

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

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran)
add_custom_command(
    OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/fortran/xcfun_autogen.F90
    DEPENDS
    genfortran
    COMMAND
    genfortran
    )
add_executable(
    genfortran
    ${CMAKE_CURRENT_SOURCE_DIR}/src/genfortran.cpp
    )
target_link_libraries(
    genfortran
    xcfun
    )
set(GENERATED_FILES
    ${CMAKE_CURRENT_BINARY_DIR}/fortran/xcfun_autogen.F90
    )
add_custom_target(xcfun_interfaces
    ALL
    DEPENDS
    ${GENERATED_FILES}
)

add_library(
    xcfun_f90_bindings
    ${CMAKE_CURRENT_BINARY_DIR}/fortran/xcfun_autogen.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran/xcfun_module.F90
    )
# Set .so version
if(BUILD_SHARED_LIBS)
 set_target_properties(xcfun_f90_bindings PROPERTIES
  VERSION ${VERSION} SOVERSION ${SOVERSION} )
endif(BUILD_SHARED_LIBS)
install(TARGETS xcfun_f90_bindings DESTINATION lib${LIB_SUFFIX})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/xcfun.mod ${CMAKE_CURRENT_BINARY_DIR}/xcfun_autogen.mod DESTINATION include)

add_executable(
    kernel_example
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran/kernel_example.F90
    )
target_link_libraries(
    kernel_example
    xcfun_f90_bindings
    xcfun
    )

include(Tests)
include(CTest)
enable_testing()
