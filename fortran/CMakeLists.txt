enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

include(FortranFlags)
include(int64)

add_library(xcfun_f90_bindings
  xcfun_module.F90
  )
set_target_properties(xcfun_f90_bindings
  PROPERTIES
    POSITION_INDEPEDENT_CODE 1
  )
target_compile_options(xcfun_f90_bindings
  PUBLIC
    "${XCFun_Fortran_FLAGS}"
    "$<$<BOOL:${OPENMP_FOUND}>:${OpenMP_Fortran_FLAGS}>"
    "$<$<BOOL:${ENABLE_64BIT_INTEGERS}>:${XCFun_64BIT_INTEGERS_FLAGS}>"
  PRIVATE
    "$<$<CONFIG:Debug>:${XCFun_Fortran_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_Fortran_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_Fortran_FLAGS_COVERAGE}>"
  )
target_compile_definitions(xcfun_f90_bindings
  PUBLIC
  "$<$<BOOL:${ENABLE_64BIT_INTEGERS}>:XCFUN_INTEGER=integer*8>"
  )

# Set .so version
if(BUILD_SHARED_LIBS)
  set_target_properties(xcfun_f90_bindings PROPERTIES VERSION
    ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION})
endif()
install(TARGETS xcfun_f90_bindings DESTINATION lib${LIB_SUFFIX})
install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/xcfun.mod DESTINATION include)

add_executable(kernel_example
  kernel_example.f90
  )
target_link_libraries(kernel_example
  xcfun_f90_bindings
  xcfun
  "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_Fortran_FLAGS_COVERAGE}>"
  )
target_compile_options(kernel_example
  PRIVATE
    "$<$<CONFIG:Debug>:${XCFun_Fortran_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_Fortran_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_Fortran_FLAGS_COVERAGE}>"
  )
set_property(TARGET kernel_example PROPERTY LINKER_LANGUAGE Fortran)
