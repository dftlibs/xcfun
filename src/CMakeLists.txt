add_library(xcfun-objlib OBJECT "")

set_target_properties(xcfun-objlib
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
target_compile_options(xcfun-objlib
  # FIXME These should be PUBLIC, but then Fortran inherits flags it doesn't
  # understand...
  PRIVATE
    "${XCFun_CXX_FLAGS}"
    "$<$<BOOL:${OPENMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
  PRIVATE
    "$<$<CONFIG:Debug>:${XCFun_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_CXX_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_CXX_FLAGS_COVERAGE}>"
  )
target_compile_definitions(xcfun-objlib
  PUBLIC
    XCFun_EXPORTS
    "$<$<BOOL:${ENABLE_64BIT_INTEGERS}>:XCFun_FORTRAN_INT=long long int>"
    XC_MAX_ORDER=${XCFun_XC_MAX_ORDER}
  )
target_include_directories(xcfun-objlib
  PUBLIC
    ${PROJECT_SOURCE_DIR}/api
    ${PROJECT_BINARY_DIR}/include
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_include_directories(xcfun-objlib
  SYSTEM
  PRIVATE
    ${PROJECT_SOURCE_DIR}/external/upstream/taylor
  )

target_sources(xcfun-objlib
  PRIVATE
    fortran.c
    xcfun.cpp
    xcint.cpp
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/config.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/densvars.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/functional.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/specmath.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xcint.hpp
  )

add_subdirectory(functionals)

if(NOT STATIC_LIBRARY_ONLY)
  add_library(xcfun-shared SHARED $<TARGET_OBJECTS:xcfun-objlib>)
  target_compile_definitions(xcfun-shared
    PUBLIC
      XCFun_EXPORTS
    INTERFACE
      $<INSTALL_INTERFACE:USING_XCFun>
    )
  target_include_directories(xcfun-shared
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/api>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  target_link_libraries(xcfun-shared
    PUBLIC
      "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:gcov>"
    )
  set_target_properties(xcfun-shared
    PROPERTIES
      VERSION ${PROJECT_VERSION_MAJOR}
      SOVERSION ${PROJECT_VERSION_MAJOR}
      CXX_VISIBILITY_PRESET hidden
      VISIBILITY_INLINES_HIDDEN 1
      MACOSX_RPATH ON
      OUTPUT_NAME "xcfun"
      EXPORT_NAME "xcfun"
    )
  install(
    TARGETS
      xcfun-shared
    EXPORT
      "XCFunTargets-shared"
    RUNTIME DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    )
endif()

if(NOT SHARED_LIBRARY_ONLY)
  add_library(xcfun-static STATIC $<TARGET_OBJECTS:xcfun-objlib>)
  target_compile_definitions(xcfun-static
    PUBLIC
      XCFun_STATIC_DEFINE
    INTERFACE
      $<INSTALL_INTERFACE:USING_XCFun>
    )
  target_include_directories(xcfun-static
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  target_link_libraries(xcfun-static
    PUBLIC
      "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:gcov>"
    )
  set_target_properties(xcfun-static
    PROPERTIES
       OUTPUT_NAME "xcfun"
       EXPORT_NAME "xcfun"
    )
  install(
    TARGETS
      xcfun-static
    EXPORT
      "XCFunTargets-static"
    ARCHIVE DESTINATION
      ${CMAKE_INSTALL_LIBDIR}
    )
endif()

if(STATIC_LIBRARY_ONLY)
  add_library(XCFun ALIAS xcfun-static)
else()
  add_library(XCFun ALIAS xcfun-shared)
endif()

include(GenerateExportHeader)
generate_export_header(XCFun
  BASE_NAME "XCFun"
  EXPORT_MACRO_NAME "XCFun_EXPORT"
  EXPORT_FILE_NAME "XCFunExport.h"
  DEPRECATED_MACRO_NAME "XCFun_DEPRECATED"
  NO_EXPORT_MACRO_NAME "XCFun_NO_EXPORT"
  STATIC_DEFINE "XCFun_STATIC_DEFINE"
  NO_DEPRECATED_MACRO_NAME "XCFun_NO_DEPRECATED"
  DEFINE_NO_DEPRECATED)
file(
  COPY ${CMAKE_CURRENT_BINARY_DIR}/XCFunExport.h
  DESTINATION ${PROJECT_BINARY_DIR}/include
  )
install(
  FILES ${PROJECT_BINARY_DIR}/include/XCFunExport.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/XCFun
  )

# <<<  Export Config  >>>
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/XCFunConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/XCFunConfig.cmake.in
  ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/XCFunConfig.cmake
  INSTALL_DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/XCFunConfig.cmake
    ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/XCFunConfigVersion.cmake
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )

if(NOT STATIC_LIBRARY_ONLY)
  install(
    EXPORT
      "XCFunTargets-shared"
    NAMESPACE
      "XCFun::"
    DESTINATION
      ${CMAKECONFIG_INSTALL_DIR}
    )
endif()

if(NOT SHARED_LIBRARY_ONLY)
  install(
    EXPORT
      "XCFunTargets-static"
    NAMESPACE
      "XCFun::"
    DESTINATION
      ${CMAKECONFIG_INSTALL_DIR}
    )
endif()

install(
  FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/config.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/densvars.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/functional.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/specmath.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xcint.hpp
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/XCFun
  )
