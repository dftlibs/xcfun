cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(XCFun CXX C Fortran)

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
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/taylor
    ${PROJECT_SOURCE_DIR}/src/functionals
    )

set(BASE_SOURCES
    src/xcint.cpp
    src/xcfun.cpp
    src/fortran.c
    )

file(GLOB FUN_SOURCES src/functionals/*.cpp)

add_library(
    xcfun
    ${BASE_SOURCES} 
    ${FUN_SOURCES} 
    )





add_executable(
    testall
    test/testall.c
    )
target_link_libraries(
    testall
    xcfun
    )


#todo:

#genfortran: src/genfortran.cpp lib/libxcfun.a
#        $(CXX) $(CXXFLAGS) $^ -o $@ -lxcfun -lm $(LIBS)
#
#fortran/xcfun_autogen.F90: genfortran include/xcfun.h
#        ./genfortran
#
#fortran/xcfun_autogen.o: fortran/xcfun_autogen.F90
#        $(FC) $(FFLAGS) -c -o $@ $^
#
#fortran/xcfun_module.o: fortran/xcfun_module.F90 fortran/xcfun_autogen.o fortran/xcfun_autogen.F90
#        $(FC) $(FFLAGS) -c -o $@ $<
#
#example_fortran: fortran/xcfun_autogen.o fortran/xcfun_module.o fortran/example.o lib/libxcfun.a
#        $(FC) $(FFLAGS) $(LIBS) -o $@ $^ -lxcfun -lm
#
#kernel_example: fortran/xcfun_autogen.o fortran/kernel_example.o fortran/xcfun_module.o  lib/libxcfun.a
#        $(FC) $(FFLAGS) $(LIBS) -o $@ $^ -lxcfun -lm


