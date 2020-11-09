add_library(xcfun
  XCFunctional.cpp
  xcint.cpp
  )

add_subdirectory(functionals)

include(GenerateExportHeader)
generate_export_header(xcfun
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/XCFunExport.h"
  DEFINE_NO_DEPRECATED
  )

target_compile_options(xcfun
  PRIVATE
    "${XCFun_CXX_FLAGS}"
    "$<$<CONFIG:Debug>:${XCFun_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_CXX_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_CXX_FLAGS_COVERAGE}>"
  )

target_compile_definitions(xcfun
  PRIVATE
    XCFUN_MAX_ORDER=${XCFUN_MAX_ORDER}
  INTERFACE
    $<INSTALL_INTERFACE:USING_XCFun>
  PUBLIC
    $<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>
  )

target_include_directories(xcfun
  PRIVATE
    ${PROJECT_SOURCE_DIR}/api
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

target_include_directories(xcfun
  SYSTEM
  PRIVATE
    ${PROJECT_SOURCE_DIR}/external/upstream/taylor
  )

target_link_libraries(xcfun
  PUBLIC
    "$<BUILD_INTERFACE:$<$<BOOL:${ENABLE_CODE_COVERAGE}>:gcov>>"
  INTERFACE
     $<INSTALL_INTERFACE:${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}>
  )

set_property(GLOBAL APPEND
  PROPERTY public_headers_property
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/XCFunExport.h
    ${CMAKE_CURRENT_SOURCE_DIR}/config.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/densvars.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/functional.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/specmath.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xcint.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/XCFunctional.hpp
  )

get_property(public_headers GLOBAL PROPERTY public_headers_property)

set_target_properties(xcfun
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
    VERSION ${PROJECT_VERSION_MAJOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    MACOSX_RPATH ON
    OUTPUT_NAME "xcfun"
    PUBLIC_HEADER "${public_headers}"
  )

install(
  TARGETS
    xcfun
  EXPORT
    "XCFunTargets"
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  PUBLIC_HEADER
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/XCFun
    COMPONENT dev
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )

# <<<  Export Config  >>>
include(CMakePackageConfigHelpers)

if(NOT DEFINED DEF_INSTALL_CMAKEDIR)
  set(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
endif()
set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/${INSTALL_CMAKEDIR}/XCFunConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/XCFunConfig.cmake.in
  ${PROJECT_BINARY_DIR}/${INSTALL_CMAKEDIR}/XCFunConfig.cmake
  INSTALL_DESTINATION
    ${INSTALL_CMAKEDIR}
  )

install(
  FILES
    ${PROJECT_BINARY_DIR}/${INSTALL_CMAKEDIR}/XCFunConfig.cmake
    ${PROJECT_BINARY_DIR}/${INSTALL_CMAKEDIR}/XCFunConfigVersion.cmake
  DESTINATION
    ${INSTALL_CMAKEDIR}
  )

install(
  EXPORT
    "XCFunTargets"
  NAMESPACE
    "XCFun::"
  DESTINATION
    ${INSTALL_CMAKEDIR}
  )

# The alias is useful with FetchContent
add_library(XCFun::xcfun ALIAS xcfun)
