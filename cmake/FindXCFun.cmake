# - Find the XCFun includes and library
#
# This module defines
#  XCFUN_INCLUDE_DIRS, where to find xcfun.h, etc.
#  XCFUN_LIBRARIES, the libraries to link against to use XCFun.
#  XCFUN_DEFINITIONS - You should add_definitons(${XCFUN_DEFINITIONS}) before
#  compiling
#  XCFUN_FOUND, If false, do not try to use XCfun.
# also defined, but not for general use are
# None of the above will be defined unles XCFun can be found.
# 

#=============================================================================
# Copyright 2010 Jonas Juselius <jonas.juselius@uit.no>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

if (XCFUN_INCLUDE_DIRS AND XCFUN_LIBRARIES)
	set(XCFUN_FIND_QUIETLY TRUE)
endif ()

find_path(XCFUN_INCLUDE_DIRS
  NAMES XCFun.h
  PATHS ${XCFUN_ROOT_DIR} $ENV{XCFUN_ROOT_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
)
find_path(XCFUN_INCLUDE_DIRS XCFun.h)

find_path(XCFUN_LIBRARIES xcfun
  PATHS ${XCFUN_ROOT_DIR} $ENV{XCFUN_ROOT_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
)
find_library(XCFUN_LIBRARIES xcfun)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XCFUN DEFAULT_MSG
                                  XCFUN_INCLUDE_DIR XCFUN_LIBRARIES)

mark_as_advanced(XCFUN_INCLUDE_DIR XCFUN_LIBRARIES)
