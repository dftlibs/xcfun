if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    set(CMAKE_Fortran_FLAGS         "-g")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "-g -traceback")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
endif()
