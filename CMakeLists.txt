cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Debug)

project(XCFun CXX C Fortran)

set (XCFun_VERSION_MAJOR 1)
set (XCFun_VERSION_MINOR 0)

option(XCFUN_NO_STDC++ "Disable libsrdc++" OFF)
option(XCFUN_ENABLE_TESTS "Disable libsrdc++" ON)

include(ConfigDefaultCompilerFlags)
include_directories(${PROJECT_BINARY_DIR})

add_subdirectory(src)
if (XCFUN_ENABLE_TESTS)
	enable_testing()
	add_subdirectory(fortran)
	add_subdirectory(test)
endif()

find_package(Doxygen)
if (DOXYGEN_FOUND)
	add_custom_target(docs
		COMMAND ${DOXYGEN_EXECUTABLE}
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		)
endif()

configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h"
	)

install(FILES include/xcfun.h DESTINATION include)

include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
	"XCFun")
set(CPACK_PACKAGE_VENDOR "Ulf Ekström")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/INSTALL")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
	"XCFun ${MRCPP_VERSION_MAJOR}.${MRCPP_VERSION_MINOR}")
include(CPack)
