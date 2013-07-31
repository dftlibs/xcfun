cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(XCFun CXX C Fortran)

# XCFun version
set(VERSION 2.0.0)
# Library version
set(SOVERSION 2.0.0)

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers" OFF)

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

include(FortranCompilers)
include(CCompilers)
include(CXXCompilers)

add_library(
    xcfun_f90_bindings
    ${CMAKE_CURRENT_SOURCE_DIR}/fortran/xcfun_module.F90
    )

# Set .so version
if(BUILD_SHARED_LIBS)
 set_target_properties(xcfun_f90_bindings PROPERTIES
  VERSION ${VERSION} SOVERSION ${SOVERSION} )
endif(BUILD_SHARED_LIBS)
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

include(Tests)
include(CTest)
enable_testing()
