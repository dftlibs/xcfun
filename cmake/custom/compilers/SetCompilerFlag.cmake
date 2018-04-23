# Adapted from
# https://github.com/robertodr/ddPCM/blob/expose-C-api/cmake/custom/compilers/SetCompilerFlag.cmake
# which was adapted by Roberto Di Remigio from
# https://github.com/SethMMorton/cmake_fortran_template/blob/master/cmake/Modules/SetCompileFlag.cmake

# Given a list of flags, this stateless function will try each, one at a time,
# and set result to the first flag that works.
# If none of the flags works, result is "".
# If the REQUIRED key is given and no flag is found, a FATAL_ERROR is raised.
#
# Call is:
# set_compile_flag(result (Fortran|C|CXX) <REQUIRED> flag1 flag2 ...)
#
# Example:
# set_compiler_flag(working_compile_flag C REQUIRED "-Wall" "-warn all")

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CheckFortranCompilerFlag)
include(CMakeParseArguments)

function(set_compiler_flag)
  set(options REQUIRED)
  set(oneValueArgs RESULT LANGUAGE)
  set(multiValueArgs FLAGS)
  cmake_parse_arguments(set_compiler_flag
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
  )

  # Silently check compiler flags
  set(restore_CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET})
  set(CMAKE_REQUIRED_QUIET TRUE)
  set(_flag_found FALSE)
  # loop over all flags, try to find the first which works
  foreach(flag IN LISTS set_compiler_flag_FLAGS)
    unset(_flag_works CACHE)
    if(${set_compiler_flag_LANGUAGE} STREQUAL "C")
      check_c_compiler_flag("${flag}" _flag_works)
    elseif(${set_compiler_flag_LANGUAGE} STREQUAL "CXX")
      check_cxx_compiler_flag("${flag}" _flag_works)
    elseif(${set_compiler_flag_LANGUAGE} STREQUAL "Fortran")
      check_Fortran_compiler_flag("${flag}" _flag_works)
    else()
      message(FATAL_ERROR "Unknown language in set_compiler_flag: ${set_compiler_flag_LANGUAGE}")
    endif()

    # if the flag works, use it, and exit loop
    # otherwise try next flag
    if(_flag_works)
      set(${set_compiler_flag_RESULT} "${flag}" PARENT_SCOPE)
      set(_flag_found TRUE)
      break()
    endif()
  endforeach()

  # raise an error if no flag was found
  if(${set_compiler_flag_REQUIRED} AND NOT _flag_found)
    message(FATAL_ERROR "None of the required flags were supported")
  endif()
  # Restore CMAKE_REQUIRED_QUIET
  set(CMAKE_REQUIRED_QUIET ${restore_CMAKE_REQUIRED_QUIET})
endfunction()
