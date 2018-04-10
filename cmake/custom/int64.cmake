# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE
# Simplified and adapted by Roberto Di Remigio

option_with_print(
  NAME
    ENABLE_64BIT_INTEGERS
  MESSAGE
    "Enable 64-bit integers"
  DEFAULT
    OFF
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
  set(XCFun_INTEGER "integer(8)")
else()
  set(XCFun_INTEGER "integer(4)")
endif()
