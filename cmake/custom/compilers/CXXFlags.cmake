# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(XCFun_CXX_FLAGS)
set(XCFun_CXX_FLAGS_DEBUG)
set(XCFun_CXX_FLAGS_RELEASE)

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
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  list(APPEND XCFun_CXX_FLAGS
    "-fno-rtti"
    "-fno-exceptions"
    "-m64"
    "-Qunused-arguments"
    "-fcolor-diagnostics"
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
    "-vec-report0"
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
    "-c9x"
    )
  list(APPEND XCFun_CXX_FLAGS_RELEASE
    "-O3"
    "-fast"
    "-Munroll"
    "-Mvect=idiom"
    "-c9x"
    )
endif()

if(ENABLE_CODE_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?[Cc]lang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3)
      message(FATAL_ERROR "Code coverage analysis on Mac OS X requires Clang version 3.0.0 or greater!")
    else()
      list(APPEND XCFun_CXX_FLAGS_DEBUG
        "-fprofile-arcs"
        "-ftest-coverage"
        )
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    list(APPEND XCFun_CXX_FLAGS_DEBUG
      "-fprofile-arcs"
      "-ftest-coverage"
      )
  else()
    message(FATAL_ERROR "Code coverage analysis requires the GNU C++ compiler!")
  endif()
endif()
