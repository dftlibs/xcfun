configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/xcfun.h
  ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/xcfun.h
  COPYONLY
  )

if(NOT DEFINED XCFun_Fortran_SOURCES)
  list(APPEND XCFun_Fortran_SOURCES
    ${PROJECT_SOURCE_DIR}/api/xcfun.f90
    )
  set(XCFun_Fortran_SOURCES "${XCFun_Fortran_SOURCES}" CACHE INTERNAL "XCFun Fortran sources for hosts")
endif()

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/xcfun.f90
  ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/xcfun.f90
  COPYONLY
  )

set_property(GLOBAL APPEND
  PROPERTY public_headers_property
    ${CMAKE_CURRENT_LIST_DIR}/xcfun.h
    ${XCFun_Fortran_SOURCES}
  )
