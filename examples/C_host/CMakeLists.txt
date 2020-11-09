cmake_minimum_required(VERSION 3.14)

project(C_host LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

find_package(XCFun CONFIG QUIET)
if(TARGET XCFun::xcfun)
  get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
  message(STATUS "Found XCFun: ${_loc} (found version ${XCFun_VERSION})")
else()
  message(STATUS "Suitable XCFun could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(xcfun
    QUIET
    URL
      https://github.com/dftlibs/xcfun/archive/v2.1.0.tar.gz
    )

  set(ENABLE_TESTALL FALSE CACHE BOOL "")
  set(XCFUN_PYTHON_INTERFACE FALSE CACHE BOOL "")
  FetchContent_MakeAvailable(xcfun)
endif()

enable_testing()

foreach(_c_src example funeval)
  add_executable(${_c_src} ${_c_src}.c)
  target_link_libraries(${_c_src}
    PUBLIC
      XCFun::xcfun
    )
  add_test(
    NAME ${_c_src}
    COMMAND $<TARGET_FILE:${_c_src}>
    )
endforeach()
