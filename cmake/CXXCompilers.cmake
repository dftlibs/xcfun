if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "-ffloat-store -fno-rtti -fno-exceptions")
    if(DEVELOPMENT_CODE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    else()
        # suppress warnings in exported code
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        set(CMAKE_CXX_FLAGS
#            "${CMAKE_CXX_FLAGS} -m32"
            "${CMAKE_CXX_FLAGS} -m64"
            )
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -m64"
            )
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize -Wno-unused")
    set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g -pg")
    if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        # radovan: vpotdamp code needs this
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
    if (ENABLE_CODE_COVERAGE)
        set (CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
        set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
    endif()
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -fopenmp"
            )
    endif()
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -static -fpic"
            )
    endif()
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177 -fno-rtti -fno-exceptions")
    if(DEVELOPMENT_CODE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    else()
        # suppress warnings in exported code
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ip")
    set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g -pg")
    set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")

    if(DEFINED MKL_FLAG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_FLAG}")
    endif()

    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -openmp"
            )
    endif()
endif ()

if(CMAKE_CXX_COMPILER_ID MATCHES PGI)
    set(CMAKE_CXX_FLAGS         "-Mpreprocess --diag_suppress 236 --diag_suppress 175") #236 suppress assert warnings and 175 suppress subscript out of range warning /SR
    set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0 -c9x")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fast -Munroll -Mvect=idiom -c9x")
    set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -g -pg")
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -mp"
            )
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES XL)
    set(CMAKE_CXX_FLAGS         " ")
    set(CMAKE_CXX_FLAGS_DEBUG   " ")
    set(CMAKE_CXX_FLAGS_RELEASE " ")
    set(CMAKE_CXX_FLAGS_PROFILE " ")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Cray)
    set(CMAKE_CXX_FLAGS         "-eZ")
    set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_CXX_FLAGS_RELEASE " ")
    set(CMAKE_CXX_FLAGS_PROFILE "-g")
endif()

if(DEFINED EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif()

