# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE
# Simplified and adapted by Roberto Di Remigio

#.rst:
#
# Enables 64-bit integer support for Fortran projects.
#
# Variables modified (provided the corresponding language is enabled)::
#
#   CMAKE_Fortran_FLAGS
#
# autocmake.yml configuration::
#
#   docopt: "--int64 Enable 64bit integers [default: False]."
#   define: "'-DENABLE_64BIT_INTEGERS={0}'.format(arguments['--int64'])"

option_with_print(
  NAME
    ENABLE_64BIT_INTEGERS
  MESSAGE
    "Enable 64-bit integers"
  DEFAULT OFF
  )

if(ENABLE_64BIT_INTEGERS)
  if(DEFINED CMAKE_Fortran_COMPILER_ID)
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
      set(XCFun_64BIT_INTEGERS_FLAGS "-fdefault-integer-8")
    endif()
    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
      set(XCFun_64BIT_INTEGERS_FLAGS "-i8")
    endif()
    if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
      set(XCFun_64BIT_INTEGERS_FLAGS "-i8")
    endif()
    if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
      set(XCFun_64BIT_INTEGERS_FLAGS "-qintsize=8 -q64")
    endif()
  endif()
endif()
