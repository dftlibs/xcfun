find_package(SWIG REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

set_property(TARGET xcfun PROPERTY POSITION_INDEPENDENT_CODE ON)

execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_SOURCE_DIR}/python")

set(XCFUN_SWIG         "${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun_swig.i")
set(XCFUN_SWIG_WRAP    "${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swigPYTHON_wrap.c")
set(CMAKE_SWIG_OUTDIR   ${CMAKE_CURRENT_BINARY_DIR}/python)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swigPYTHON_wrap.c
                   COMMAND "${SWIG_EXECUTABLE}"
                   ARGS "-python" "-I${CMAKE_CURRENT_SOURCE_DIR}/include" "-outdir" "${CMAKE_SWIG_OUTDIR}" "-o" ${XCFUN_SWIG_WRAP} ${XCFUN_SWIG}
                   DEPENDS xcfun ${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun_swig.i)

set(XCFUN_PY       "${CMAKE_CURRENT_SOURCE_DIR}/python/xcfun.py")
set(XCFUN_SWIG_PY  "${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swig.py")
set(XCFUN_PY_LIB   "${CMAKE_CURRENT_BINARY_DIR}/python/_xcfun_swig.so")

set(SETUP_PY_IN  "${CMAKE_CURRENT_SOURCE_DIR}/python/setup.py.in")
set(SETUP_PY     "${CMAKE_CURRENT_BINARY_DIR}/python/setup.py")

configure_file(${SETUP_PY_IN} ${SETUP_PY})

add_custom_command(OUTPUT ${XCFUN_PY_LIB}
                   COMMAND ${PYTHON_EXECUTABLE}
                   ARGS "setup.py" "build_ext" "--inplace"
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python
                   DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/python/setup.py ${CMAKE_CURRENT_BINARY_DIR}/python/xcfun_swigPYTHON_wrap.c)

add_custom_target(target ALL DEPENDS ${XCFUN_PY_LIB})

install(FILES "${XCFUN_PY}"           DESTINATION "python")
install(FILES "${XCFUN_SWIG_PY}"      DESTINATION "python")
install(FILES "${XCFUN_PY_LIB}"       DESTINATION "python")
