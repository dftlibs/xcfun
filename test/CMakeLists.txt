add_executable(testall
  testall.c
  )
target_link_libraries(testall
  XCFun
  )
target_compile_options(testall
  PRIVATE
    "${XCFun_C_FLAGS}"
    "$<$<BOOL:${OPENMP_FOUND}>:${OpenMP_C_FLAGS}>"
    "$<$<CONFIG:Debug>:${XCFun_C_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_C_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_C_FLAGS_COVERAGE}>"
  )
add_test(testall ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testall)

if(TEST_Fortran_BINDINGS)
  add_executable(kernel_example
    kernel_example.f90
    ${PROJECT_SOURCE_DIR}/api/xcfun.F90
    )
  target_link_libraries(kernel_example
    XCFun
    )
  target_compile_options(kernel_example
    PRIVATE
      "$<$<CONFIG:Debug>:${XCFun_Fortran_FLAGS_DEBUG}>"
      "$<$<CONFIG:Release>:${XCFun_Fortran_FLAGS_RELEASE}>"
      "$<$<BOOL:${ENABLE_64BIT_INTEGERS}>:${XCFun_64BIT_INTEGERS_FLAGS}>"
      "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_Fortran_FLAGS_COVERAGE}>"
    )
  set_property(TARGET kernel_example PROPERTY LINKER_LANGUAGE Fortran)
  add_test(kernel_example ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kernel_example)
endif()

if(ENABLE_PYTHON_INTERFACE)
  add_test(
    NAME python-interface
    COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/test/test.py"
    ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH}
    )
endif()

