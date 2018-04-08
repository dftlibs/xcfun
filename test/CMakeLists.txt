add_executable(testall
  testall.c
  )
target_link_libraries(
  testall
  xcfun
  )
target_compile_options(testall
  PRIVATE
    "${XCFun_C_FLAGS}"
    "$<$<BOOL:${OPENMP_FOUND}>:${OpenMP_C_FLAGS}>"
    "$<$<CONFIG:Debug>:${XCFun_C_FLAGS_DEBUG}>"
    "$<$<CONFIG:Debug>:${XCFun_C_FLAGS_DEBUG}>"
  )
add_test(testall ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testall)

if(ENABLE_Fortran_INTERFACE)
  add_test(kernel_example ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kernel_example)
endif()

if(ENABLE_PYTHON_INTERFACE)
  add_test(
    NAME python-interface
    COMMAND "${PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/test/test.py"
    ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH}
    )
endif()
