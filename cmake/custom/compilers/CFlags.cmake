# C99
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(XCFun_C_FLAGS)
set(XCFun_C_FLAGS_DEBUG)
set(XCFun_C_FLAGS_RELEASE)

if(CMAKE_C_COMPILER_ID MATCHES GNU)
  list(APPEND XCFun_C_FLAGS
    "-ffloat-store"
    "-m64"
    )
  list(APPEND XCFun_C_FLAGS_DEBUG
    "-Wall"
    "-O0"
    "-g3"
    "-Wextra"
    "-Winit-self"
    "-Wuninitialized"
    "-Wmissing-declarations"
    "-Wwrite-strings"
    "-Wno-sign-compare"
    )
  list(APPEND XCFun_C_FLAGS_RELEASE
    "-O3"
    "-ffast-math"
    "-funroll-loops"
    "-ftree-vectorize"
    "-Wno-unused"
    )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  list(APPEND XCFun_C_FLAGS
    "-m64"
    "-Qunused-arguments"
    "-fcolor-diagnostics"
    )
  list(APPEND XCFun_C_FLAGS_DEBUG
    "-Wall"
    "-O0"
    "-g3"
    "-Wextra"
    "-Winit-self"
    "-Wuninitialized"
    "-Wmissing-declarations"
    "-Wwrite-strings"
    "-Wno-sign-compare"
    )
  list(APPEND XCFun_C_FLAGS_RELEASE
    "-O3"
    "-ffast-math"
    "-funroll-loops"
    "-ftree-vectorize"
    "-Wno-unused"
    )
endif()

if(CMAKE_C_COMPILER_ID MATCHES Intel)
  list(APPEND XCFun_C_FLAGS
    "-g"
    "-wd981"
    "-wd279"
    "-wd383"
    "-vec-report0"
    "-wd1572"
    "-wd1777"
    "-restrict"
    )
  list(APPEND XCFun_C_FLAGS_DEBUG
    "-Wall"
    "-O0"
    "-g"
    "-w3"
    "-Wuninitialized"
    "-Wno-sign-compare"
    )
  list(APPEND XCFun_C_FLAGS_RELEASE
    "-O3"
    "-ip"
    )
endif ()

if(CMAKE_C_COMPILER_ID MATCHES PGI)
  list(APPEND XCFun_C_FLAGS
    "-Mpreprocess"
    )
  list(APPEND XCFun_C_FLAGS_DEBUG
    "-g"
    "-O0"
    "-c9x"
    )
  list(APPEND XCFun_C_FLAGS_RELEASE
    "-O3"
    "-fast"
    "-Munroll"
    "-Mvect=idiom"
    "-c9x"
    )
endif()