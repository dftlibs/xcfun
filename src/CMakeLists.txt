add_library(xcfun OBJECT "")

include(GenerateExportHeader)
generate_export_header(xcfun
  BASE_NAME "XCFun"
  EXPORT_MACRO_NAME "XCFun_EXPORT"
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/XCFunExport.h"
  DEPRECATED_MACRO_NAME "XCFun_DEPRECATED"
  NO_EXPORT_MACRO_NAME "XCFun_NO_EXPORT"
  STATIC_DEFINE "XCFun_STATIC_DEFINE"
  NO_DEPRECATED_MACRO_NAME "XCFun_NO_DEPRECATED"
  DEFINE_NO_DEPRECATED
  )
install(
  FILES ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/XCFunExport.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/XCFun
  )

set_target_properties(xcfun
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )
target_compile_options(xcfun
  # FIXME These should be PUBLIC, but then Fortran inherits flags it doesn't
  # understand...
  PRIVATE
    "${XCFun_CXX_FLAGS}"
  PRIVATE
    "$<$<CONFIG:Debug>:${XCFun_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_CXX_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_CXX_FLAGS_COVERAGE}>"
  )
target_compile_definitions(xcfun
  PUBLIC
    XCFun_EXPORTS
    XC_MAX_ORDER=${XCFun_XC_MAX_ORDER}
  )
target_include_directories(xcfun
  PRIVATE
    ${PROJECT_SOURCE_DIR}/api
    ${PROJECT_BINARY_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_include_directories(xcfun
  SYSTEM
  PRIVATE
    ${PROJECT_SOURCE_DIR}/external/upstream/taylor
  )

target_sources(xcfun
  PRIVATE
    xcfun.cpp
    xcint.cpp
  )

add_library(xcfun-c-objlib OBJECT fortran.c)
set_target_properties(xcfun-c-objlib
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
  )
target_compile_options(xcfun-c-objlib
  # FIXME These should be PUBLIC, but then Fortran inherits flags it doesn't
  # understand...
  PRIVATE
    "${XCFun_C_FLAGS}"
  PRIVATE
    "$<$<CONFIG:Debug>:${XCFun_C_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_C_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_C_FLAGS_COVERAGE}>"
  )
target_compile_definitions(xcfun-c-objlib
  PUBLIC
    XCFun_EXPORTS
    "$<$<BOOL:${ENABLE_64BIT_INTEGERS}>:XCFun_FORTRAN_INT=long long int>"
  )
target_include_directories(xcfun-c-objlib
  PRIVATE
    ${PROJECT_SOURCE_DIR}/api
    ${PROJECT_BINARY_DIR}/include
  )

add_subdirectory(functionals)

if(NOT STATIC_LIBRARY_ONLY)
  add_library(xcfun-shared SHARED
    $<TARGET_OBJECTS:xcfun>
    $<TARGET_OBJECTS:xcfun-c-objlib>
    )
  target_compile_definitions(xcfun-shared
    PUBLIC
      XCFun_EXPORTS
    INTERFACE
      $<INSTALL_INTERFACE:USING_XCFun>
    )
  target_include_directories(xcfun-shared
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
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
  add_library(xcfun-static STATIC
    $<TARGET_OBJECTS:xcfun>
    $<TARGET_OBJECTS:xcfun-c-objlib>
    )
  target_compile_definitions(xcfun-static
    PUBLIC
      XCFun_STATIC_DEFINE
    INTERFACE
      $<INSTALL_INTERFACE:USING_XCFun>
    )
  target_include_directories(xcfun-static
    INTERFACE
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
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

# The second alias is useful with FetchContent
if(STATIC_LIBRARY_ONLY)
  add_library(XCFun ALIAS xcfun-static)
  add_library(XCFun::xcfun ALIAS xcfun-static)
else()
  add_library(XCFun ALIAS xcfun-shared)
  add_library(XCFun::xcfun ALIAS xcfun-shared)
endif()

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
