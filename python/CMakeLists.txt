set(Python_ADDITIONAL_VERSIONS 3.7 3.6 3.5)  # adjust with CMake minimum FindPythonInterp
find_package(PythonLibsNew 2.7 REQUIRED)
find_package(SWIG REQUIRED)
find_package(NumPy REQUIRED)

set(XCFUN_SWIG ${CMAKE_CURRENT_LIST_DIR}/xcfun_swig.i)
set(XCFUN_SWIG_WRAP "xcfun_swigPYTHON_wrap.c")

add_custom_command(
  OUTPUT ${XCFUN_SWIG_WRAP}
  COMMAND "${SWIG_EXECUTABLE}"
  ARGS "-python" "-I${PROJECT_SOURCE_DIR}/include" "-outdir" "${CMAKE_CURRENT_BINARY_DIR}" "-o" ${XCFUN_SWIG_WRAP} ${XCFUN_SWIG}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS xcfun ${XCFUN_SWIG}
  )

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/setup.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/setup.py
  )

set(XCFUN_PY "${CMAKE_CURRENT_LIST_DIR}/xcfun.py")
set(XCFUN_SWIG_PY "${CMAKE_CURRENT_BINARY_DIR}/xcfun_swig.py")
set(XCFUN_PY_LIB "${CMAKE_CURRENT_BINARY_DIR}/_xcfun_swig.so")

add_custom_command(
  OUTPUT ${XCFUN_PY_LIB}
  COMMAND ${PYTHON_EXECUTABLE}
  ARGS "setup.py" "build_ext" "--inplace"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/setup.py ${CMAKE_CURRENT_BINARY_DIR}/xcfun_swigPYTHON_wrap.c
  )

add_custom_target(target ALL DEPENDS ${XCFUN_PY_LIB})

install(FILES "${XCFUN_PY}"           DESTINATION "python")
install(FILES "${XCFUN_SWIG_PY}"      DESTINATION "python")
install(FILES "${XCFUN_PY_LIB}"       DESTINATION "python")
