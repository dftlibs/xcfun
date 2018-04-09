set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(BASE_SOURCES
    xcint.cpp
    xcfun.cpp
    fortran.c
    )

file(GLOB FUNCTIONAL_SOURCES functionals/*.cpp)

add_library(
    xcfun
    ${BASE_SOURCES}
    ${FUNCTIONAL_SOURCES}
    )
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

install(TARGETS xcfun DESTINATION lib${LIB_SUFFIX})
file(GLOB headers "include/*.h")
install(FILES ${headers} DESTINATION include)
