# FIXME Remove when we finally become Python 3 only
set(PYBIND11_PYTHON_VERSION 2.7)

if(MSVC)
  set(PYBIND11_CPP_STANDARD "/std:c++${CMAKE_CXX_STANDARD}")
else()
  set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")
endif()

set(PYBIND11_VERSION 2.3.0)

find_package(pybind11 ${PYBIND11_VERSION} CONFIG QUIET)
if(pybind11_FOUND)
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
  message(STATUS "Suitable pybind11 could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Populate(pybind11_sources
    QUIET
    URL
      https://github.com/pybind/pybind11/archive/v${PYBIND11_VERSION}.tar.gz
    )

  set(PYBIND11_PYTHON_VERSION ${PYBIND11_PYTHON_VERSION})
  set(PYBIND11_TEST OFF CACHE BOOL "")
  set(PYMOD_INSTALL_FULLDIR ${PYMOD_INSTALL_FULLDIR})

  add_subdirectory(
    ${pybind11_sources_SOURCE_DIR}
    ${pybind11_sources_BINARY_DIR}
    )
endif()

