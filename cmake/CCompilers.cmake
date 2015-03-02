if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS         "-ffloat-store")
    if(DEVELOPMENT_CODE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
    else()
        # suppress warnings in exported code
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        set(CMAKE_C_FLAGS
#            "${CMAKE_C_FLAGS} -m32"
            "${CMAKE_C_FLAGS} -m64"
            )
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -m64"
            )
    endif()
    set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize -Wno-unused")
    set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")
    if (ENABLE_CODE_COVERAGE)
        set (CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
        set (CMAKE_C_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
    endif()
    if(ENABLE_OMP)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -fopenmp"
            )
    endif()
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -static -fpic"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd1777 -restrict")
    if(NOT DEVELOPMENT_CODE)
        # suppress warnings in exported code
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    endif()
    set(CMAKE_C_FLAGS_DEBUG   "-O0")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ip")
    set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")

    if(DEFINED MKL_FLAG)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MKL_FLAG}")
    endif()

    if(ENABLE_OMP)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -openmp"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES PGI)
    set(CMAKE_C_FLAGS         "-Mpreprocess")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0 -c9x")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -fast -Munroll -Mvect=idiom -c9x")
    set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")
    if(ENABLE_OMP)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -mp"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES XL)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   " ")
    set(CMAKE_C_FLAGS_RELEASE " ")
    set(CMAKE_C_FLAGS_PROFILE " ")
endif()

if(CMAKE_C_COMPILER_ID MATCHES Cray)
    set(CMAKE_C_FLAGS         "-eZ")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE " ")
    set(CMAKE_C_FLAGS_PROFILE "-g")
endif()

if(DEFINED EXTRA_C_FLAGS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_C_FLAGS}")
endif()

