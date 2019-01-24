find_package(SWIG REQUIRED)
include(UseSWIG)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/custom)
find_package(PythonLibsNew 2.7 REQUIRED)
find_package(NumPy REQUIRED)

set_source_files_properties(xcfun_swig.i
  PROPERTIES
    CPLUSPLUS
      ON
    SWIG_FLAGS
      "-I${PROJECT_BINARY_DIR}/include;-I${PROJECT_SOURCE_DIR}/api"
  )
set(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})

if(CMAKE_VERSION VERSION_LESS "3.8.0")
  swig_add_module(xcfun_swig python xcfun_swig.i)
else()
  swig_add_library(xcfun_swig LANGUAGE python SOURCES xcfun_swig.i)
endif()

target_link_libraries(${SWIG_MODULE_xcfun_swig_REAL_NAME}
  PRIVATE
    ${PYTHON_LIBRARIES}
    XCFun
  )
target_include_directories(${SWIG_MODULE_xcfun_swig_REAL_NAME}
  PRIVATE
    ${PYTHON_INCLUDE_DIRS}
    ${NUMPY_INCLUDE_DIRS}
  )
file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR}
                        ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(XCFun_RPATH "@loader_path/${_rel}")
else()
  set(XCFun_RPATH "\$ORIGIN/${_rel}")
endif()
set_target_properties(${SWIG_MODULE_xcfun_swig_REAL_NAME}
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${XCFun_RPATH}${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH ON
  )
set_target_properties(${SWIG_MODULE_xcfun_swig_REAL_NAME}
  PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  )

configure_file(init.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py COPYONLY)
file(
  COPY
    xcfun.py
  DESTINATION
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/xcfun.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/xcfun_swig.py
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )
install(
  TARGETS
    ${SWIG_MODULE_xcfun_swig_REAL_NAME}
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )
