add_executable(testall
  testall.c
  )
target_link_libraries(
  testall
  xcfun
  )
# FIXME
target_include_directories(testall
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
  )
add_test(testall ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testall)

if(ENABLE_Fortran_INTERFACE)
  add_test(kernel_example ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kernel_example)
endif()

if(ENABLE_PYTHON_INTERFACE)
  add_test(NAME python-interface  COMMAND "${PYTHON_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/test/test.py")
endif()
