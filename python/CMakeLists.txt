# create Python module
pybind11_add_module(_xcfun MODULE THIN_LTO export_xcfun.cpp)

target_link_libraries(_xcfun
  PRIVATE
    XCFun
  )

target_include_directories(_xcfun
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/external/upstream/taylor
  )

set_target_properties(_xcfun
  PROPERTIES
  PREFIX "${PYTHON_MODULE_PREFIX}"
  SUFFIX "${PYTHON_MODULE_EXTENSION}"
  )

file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(XCFun_RPATH "@loader_path/${_rel}")
else()
  set(XCFun_RPATH "\$ORIGIN/${_rel}")
endif()
set_target_properties(_xcfun
  PROPERTIES
  MACOSX_RPATH ON
  SKIP_BUILD_RPATH OFF
  BUILD_WITH_INSTALL_RPATH OFF
  INSTALL_RPATH "${XCFun_RPATH}${CMAKE_INSTALL_LIBDIR}"
  INSTALL_RPATH_USE_LINK_PATH ON
  )
set_target_properties(_xcfun
  PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
    RESOURCE
      "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py;${CMAKE_CURRENT_SOURCE_DIR}/xcfun.py"
  )

foreach(_f __init__.py xcfun.py)
  configure_file(
    ${_f}
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/${_f}
    COPYONLY
    )
endforeach()

install(
  TARGETS
    _xcfun
  LIBRARY
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  RESOURCE
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  )
