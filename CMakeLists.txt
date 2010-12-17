cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Debug)

project(XCFun CXX C Fortran)

set (XCFun_VERSION_MAJOR 1)
set (XCFun_VERSION_MINOR 0)

option(XCFUN_NO_STDC++ "Disable libsrdc++" OFF)
option(XCFUN_ENABLE_TESTS "Enable tests" ON)

include(ConfigCompilerFlags)
if (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -fno-rtti -fno-exceptions ")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-g -O2 -DNDEBUG -Wno-unused")

	if (XCFUN_NO_STDC++)
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-threadsafe-statics")
	endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
	set (CMAKE_CXX_FLAGS "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177 -fno-rtti -fno-exceptions")
	set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
	set (CMAKE_CXX_FLAGS_RELEASE "-g -O2 -DNDEBUG")
	set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
endif ()
config_compiler_flags(C++)

if (CMAKE_COMPILER_IS_GNUC)
	set (CMAKE_C_FLAGS "-Wall")
	set (CMAKE_C_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_C_FLAGS_RELEASE "-g -O2 -DNDEBUG -Wno-unused")
elseif (CMAKE_C_COMPILER_ID MATCHES Intel)
	set (CMAKE_C_FLAGS "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
	set (CMAKE_C_FLAGS_DEBUG "-g -O0")
	set (CMAKE_C_FLAGS_RELEASE "-g -O2 -DNDEBUG")
	set (CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
endif ()
config_compiler_flags(C)

if (CMAKE_COMPILER_IS_GNUFORTRAN)
	set (CMAKE_Fortran_FLAGS "-Wall -Jfortran")
	set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
	set (CMAKE_Fortran_FLAGS_RELEASE "-g -O2")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	set (CMAKE_Fortran_FLAGS "")
	set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
	set (CMAKE_Fortran_FLAGS_RELEASE "-g -O2")
	set (CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -shared-intel")
endif ()
config_compiler_flags(Fortran)

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

find_program(VALGRIND_EXECUTABLE valgrind)
if (VALGRIND_EXECUTABLE)
	set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE})
	set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full")
	#set(MEMORYCHECK_SUPPRESSIONS_FILE ${CMAKE_BINARY_DIR}/valgrind.supp)
endif()
configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake.in
	        ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

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
