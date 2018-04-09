add_library(xcfun "")

set_target_properties(xcfun
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )
target_compile_options(xcfun
  # FIXME These should be PUBLIC, but then Fortran inherits flags it doesn't
  # understand...
  PRIVATE
    "${XCFun_CXX_FLAGS}"
    "$<$<BOOL:${OPENMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
  PRIVATE
    "$<$<CONFIG:Debug>:${XCFun_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_CXX_FLAGS_RELEASE}>"
  )
target_include_directories(xcfun
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/functionals
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_include_directories(xcfun
  SYSTEM
  PRIVATE
    ${PROJECT_SOURCE_DIR}/external/upstream/taylor
  )

# Set .so version
if(BUILD_SHARED_LIBS)
 set_target_properties(xcfun PROPERTIES
   VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION} )
endif()

target_sources(xcfun
  PRIVATE
    fortran.c
    xcfun.cpp
    xcint.cpp
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/config.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/densvars.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/functional.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/specmath.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xcint.hpp
  )

add_subdirectory(functionals)

install(TARGETS xcfun DESTINATION lib${LIB_SUFFIX})
file(GLOB headers "include/*.h")
install(FILES ${headers} DESTINATION include)
