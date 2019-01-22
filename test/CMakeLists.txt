add_executable(testall
  testall.c
  )
target_link_libraries(testall
  XCFun
  )
target_compile_options(testall
  PRIVATE
    "${XCFun_C_FLAGS}"
    "$<$<CONFIG:Debug>:${XCFun_C_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_C_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_C_FLAGS_COVERAGE}>"
  )
add_test(testall ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testall)

if(ENABLE_FC_SUPPORT)
  add_executable(fortran_example
    examples/fortran_example.f90
    ${PROJECT_SOURCE_DIR}/api/xcfun.F90
    )
  target_link_libraries(fortran_example
    XCFun
    )
  target_compile_options(fortran_example
    PRIVATE
      "$<$<CONFIG:Debug>:${XCFun_Fortran_FLAGS_DEBUG}>"
      "$<$<CONFIG:Release>:${XCFun_Fortran_FLAGS_RELEASE}>"
      "$<$<BOOL:${ENABLE_64BIT_INTEGERS}>:${XCFun_64BIT_INTEGERS_FLAGS}>"
      "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_Fortran_FLAGS_COVERAGE}>"
    )
  set_target_properties(fortran_example PROPERTIES LINKER_LANGUAGE Fortran)
  add_test(fortran_example ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/fortran_example)
endif()

if(ENABLE_PYTHON_INTERFACE)
  add_test(
    NAME python-interface
    COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/test.py
    )
  set_tests_properties(python-interface
    PROPERTIES
    ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}:$ENV{PYTHONPATH}
    )
endif()

