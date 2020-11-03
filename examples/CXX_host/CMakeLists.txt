cmake_minimum_required(VERSION 3.14)

project(CXX_host LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

foreach(_src example)
  add_executable(${_src} ${_src}.cpp)
  target_link_libraries(${_src}
    PUBLIC
      XCFun::xcfun
    )
  add_test(
    NAME ${_src}
    COMMAND $<TARGET_FILE:${_src}>
    )
endforeach()
