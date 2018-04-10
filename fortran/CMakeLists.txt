include(FortranFlags)
include(int64)

configure_file(xcfun.F90.in ${PROJECT_BINARY_DIR}/include/xcfun.F90 @ONLY)
set_source_files_properties(xcfun.F90
  PROPERTIES
    GENERATED 1
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/include/xcfun.F90
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/XCFun
  )
