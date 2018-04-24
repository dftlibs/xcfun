#.rst:
#
# Enables creation of static/shared library.
# If the shared library is created, make it as static as possible.
#
# Variables modified (provided the corresponding language is enabled)::
#
#   XCFun_Fortran_FLAGS
#   XCFun_C_FLAGS
#   XCFun_CXX_FLAGS
#
# autocmake.yml configuration::
#
#   docopt: "--static Create only the static library [default: False]."
#   define: "'-DSTATIC_LIBRARY_ONLY={0}'.format(arguments['--static'])"

option_with_print(STATIC_LIBRARY_ONLY "Create the static library only" OFF)
option_with_print(SHARED_LIBRARY_ONLY "Create the shared library only" OFF)
option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)

if(ENABLE_GENERIC)
  if(DEFINED CMAKE_Fortran_COMPILER_ID)
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
      list(APPEND XCFun_Fortran_FLAGS
        "-static-libgfortran"
        )
    endif()
    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
      list(APPEND XCFun_Fortran_FLAGS
        "-static-libgcc"
        "-static-intel"
        )
    endif()
  endif()

  if(DEFINED CMAKE_C_COMPILER_ID)
    if(CMAKE_C_COMPILER_ID MATCHES GNU)
      list(APPEND XCFun_C_FLAGS
        "-static-libgcc"
        "-fpic"
        )
    endif()
    if(CMAKE_C_COMPILER_ID MATCHES Intel)
      list(APPEND XCFun_C_FLAGS
        "-static-libgcc"
        "-static-intel"
        "-wd10237"
        )
    endif()
    if(CMAKE_C_COMPILER_ID MATCHES Clang)
      list(APPEND XCFun_C_FLAGS
        "-fpic"
        )
    endif()
  endif()

  if(DEFINED CMAKE_CXX_COMPILER_ID)
    if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
      list(APPEND XCFun_CXX_FLAGS
        "-static-libstdc++"
        "-static-libgcc"
        )
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
      list(APPEND XCFun_CXX_FLAGS
        "-Wl,--as-needed"
        "-static-libstdc++"
        "-static-libgcc"
        "-static-intel"
        "-wd10237"
        )
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
      list(APPEND XCFun_CXX_FLAGS
        "-static-libstdc++"
        )
    endif()
  endif()
endif()
