enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
add_library(xcfun_f90_bindings
  xcfun_module.F90
  )

# Set .so version
if(BUILD_SHARED_LIBS)
  set_target_properties(xcfun_f90_bindings PROPERTIES VERSION
    ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION})
endif()
install(TARGETS xcfun_f90_bindings DESTINATION lib${LIB_SUFFIX})
install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/xcfun.mod DESTINATION include)

if(ENABLE_64BIT_INTEGERS)
  add_definitions("-DXCFUN_INTEGER=integer*8")
  add_definitions("-DXCFUN_FORTRAN_INT=long long int")
endif()

add_executable(kernel_example
  kernel_example.F90
  )
target_link_libraries(
  kernel_example
  xcfun_f90_bindings
  xcfun
  )
set_property(TARGET kernel_example PROPERTY LINKER_LANGUAGE Fortran)
