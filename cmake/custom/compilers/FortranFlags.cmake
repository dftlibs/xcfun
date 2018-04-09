set(XCFun_Fortran_FLAGS)
set(XCFun_Fortran_FLAGS_DEBUG)
set(XCFun_Fortran_FLAGS_RELEASE)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  list(APPEND XCFun_Fortran_FLAGS
    "-ffloat-store"
    "-fcray-pointer"
    "-m64"
    "-fimplicit-none"
    "-fautomatic"
    "-fmax-errors=5"
    )
  list(APPEND XCFun_Fortran_FLAGS_DEBUG
    "-Wall"
    "-Wuninitialized"
    "-O0"
    "-g3"
    "-fbacktrace"
    )
  list(APPEND XCFun_Fortran_FLAGS_RELEASE
    "-O3"
    "-ffast-math"
    "-funroll-loops"
    "-ftree-vectorize"
    )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  list(APPEND XCFun_Fortran_FLAGS
    "-fpp"
    "-assume byterecl"
    "-nosave"
    )
  list(APPEND XCFun_Fortran_FLAGS_DEBUG
    "-Wall"
    "-O0"
    "-g"
    "-traceback"
    )
  list(APPEND XCFun_Fortran_FLAGS_RELEASE
    "-O3"
    "-ip"
    )
endif ()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  list(APPEND XCFun_Fortran_FLAGS
    "-mcmodel=medium"
    "-pgcpplibs"
    )
  list(APPEND XCFun_Fortran_FLAGS_DEBUG
    "-g"
    "-O0"
    "-Mframe"
    )
  list(APPEND XCFun_Fortran_FLAGS_RELEASE
    "-O3"
    "-Mipa=fast"
    "-Munroll"
    "-Mvect=idiom"
    "-c9x"
    )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  list(APPEND XCFun_Fortran_FLAGS_DEBUG
    "-fprofile-arcs"
    "-ftest-coverage"
    )
else()
  message(FATAL_ERROR "Code coverage analysis requires the GNU Fortran compiler!")
endif()
