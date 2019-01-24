set(XCFun_CXX_FLAGS)
set(XCFun_CXX_FLAGS_DEBUG)
set(XCFun_CXX_FLAGS_RELEASE)
set(XCFun_CXX_FLAGS_COVERAGE)

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  list(APPEND XCFun_CXX_FLAGS
    "-ffloat-store"
    "-fno-rtti"
    "-fno-exceptions"
    "-m64"
    )
  list(APPEND XCFun_CXX_FLAGS_DEBUG
    "-Wall"
    "-O0"
    "-g3"
    "-Wextra"
    "-Winit-self"
    "-Woverloaded-virtual"
    "-Wuninitialized"
    "-Wmissing-declarations"
    "-Wwrite-strings"
    "-Weffc++"
    "-Wno-sign-compare"
    )
  list(APPEND XCFun_CXX_FLAGS_RELEASE
    "-O3"
    "-ffast-math"
    "-funroll-loops"
    "-ftree-vectorize"
    "-Wno-unused"
    )
  list(APPEND XCFun_CXX_FLAGS_COVERAGE
    "${CODE_COVERAGE_FLAGS}"
    )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  list(APPEND XCFun_CXX_FLAGS
    "-fno-rtti"
    "-fno-exceptions"
    "-m64"
    )
  list(APPEND XCFun_CXX_FLAGS_DEBUG
    "-Wall"
    "-O0"
    "-g3"
    "-Wextra"
    "-Winit-self"
    "-Woverloaded-virtual"
    "-Wuninitialized"
    "-Wmissing-declarations"
    "-Wwrite-strings"
    "-Weffc++"
    "-Wdocumentation"
    "-Wno-sign-compare"
    )
  list(APPEND XCFun_CXX_FLAGS_RELEASE
    "-O3"
    "-ffast-math"
    "-funroll-loops"
    "-ftree-vectorize"
    "-Wno-unused"
    )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  list(APPEND XCFun_CXX_FLAGS
    "-g"
    "-wd981"
    "-wd279"
    "-wd383"
    "-wd1572"
    "-wd177"
    "-fno-rtti"
    "-fno-exceptions"
    )
  list(APPEND XCFun_CXX_FLAGS_DEBUG
    "-Wall"
    "-O0"
    )
  list(APPEND XCFun_CXX_FLAGS_RELEASE
    "-O3"
    "-ip"
    )
endif ()

if(CMAKE_CXX_COMPILER_ID MATCHES PGI)
  #236 suppress assert warnings and 175 suppress subscript out of range warning /SR
  list(APPEND XCFun_CXX_FLAGS
    "-Mpreprocess"
    "--diag_suppress 236"
    "--diag_suppress 175"
    )
  list(APPEND XCFun_CXX_FLAGS_DEBUG
    "-g"
    "-O0"
    )
  list(APPEND XCFun_CXX_FLAGS_RELEASE
    "-O3"
    "-fast"
    "-Munroll"
    "-Mvect=idiom"
    )
endif()

message(STATUS "C++ compiler flags     : ${CMAKE_CXX_FLAGS} ${XCFun_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${XCFun_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
