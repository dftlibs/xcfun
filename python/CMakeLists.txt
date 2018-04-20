find_package(SWIG REQUIRED)
INCLUDE(UseSWIG)

find_package(PythonInterp 2.7 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)
find_package(NumPy REQUIRED)

set(XCFUN_SWIG ${CMAKE_CURRENT_LIST_DIR}/xcfun_swig.i)

set_property(SOURCE "${XCFUN_SWIG}" PROPERTY SWIG_FLAGS "-I${CMAKE_BINARY_DIR}/include" 
                                                        "-I${PROJECT_SOURCE_DIR}/api")
set_property(SOURCE "${XCFUN_SWIG}" PROPERTY SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})

# requires CMake 3.8
#swig_add_library(xcfun_swig LANGUAGE python SOURCES "${XCFUN_SWIG}")
swig_add_module(xcfun_swig python "${XCFUN_SWIG}")

target_include_directories(${SWIG_MODULE_xcfun_swig_REAL_NAME} 
                           PRIVATE ${PYTHON_INCLUDE_DIRS}
                                   ${NUMPY_INCLUDE_DIRS}
                                   ${PROJECT_BINARY_DIR}/include
                                   ${PROJECT_SOURCE_DIR}/api)
target_link_libraries(${SWIG_MODULE_xcfun_swig_REAL_NAME} "${PYTHON_LIBRARIES}" xcfun-static)
set_target_properties(${SWIG_MODULE_xcfun_swig_REAL_NAME} PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(pyxcfun DEPENDS ${SWIG_MODULE_xcfun_swig_REAL_NAME})

set(XCFUN_PY "${CMAKE_CURRENT_LIST_DIR}/xcfun.py")
file(COPY "${XCFUN_PY}" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "${XCFUN_PY}" DESTINATION "python")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/xcfun_swig.py" DESTINATION "python")
install(TARGETS ${SWIG_MODULE_xcfun_swig_REAL_NAME} DESTINATION "python")
