cmake_minimum_required(VERSION 3.11)

project(Fortran_host LANGUAGES Fortran)

find_package(XCFun CONFIG QUIET)
if(TARGET XCFun::xcfun)
  get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
  message(STATUS "Found XCFun: ${_loc} (found version ${XCFun_VERSION})")
else()
  message(STATUS "Suitable XCFun could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(xcfun_sources
    QUIET
    URL
      https://github.com/dftlibs/xcfun/archive/v2.0.0a2.tar.gz
    )

  FetchContent_GetProperties(xcfun_sources)

  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external)
  set(ENABLE_TESTALL FALSE CACHE BOOL "")
  set(XCFUN_PYTHON_INTERFACE FALSE)

  if(NOT xcfun_sources_POPULATED)
    FetchContent_Populate(xcfun_sources)

    add_subdirectory(
      ${xcfun_sources_SOURCE_DIR}
      ${xcfun_sources_BINARY_DIR}
      )
  endif()
endif()

enable_testing()

foreach(_src example)
  add_executable(${_src} ${_src}.f90)
  target_sources(${_src}
    PRIVATE
      ${XCFun_Fortran_SOURCES}
    )
  target_link_libraries(${_src}
    PUBLIC
      XCFun::xcfun
    )
  add_test(
    NAME ${_src}
    COMMAND $<TARGET_FILE:${_src}>
    )
endforeach()
