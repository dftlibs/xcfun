# Compile C++ example
add_executable(cpp_example
  cpp_example.cpp
)
target_link_libraries(cpp_example
  XCFun::xcfun
)
target_compile_options(cpp_example
  PRIVATE
    "$<$<CONFIG:Debug>:${XCFun_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_CXX_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_64BIT_INTEGERS}>:${XCFun_64BIT_INTEGERS_FLAGS}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_CXX_FLAGS_COVERAGE}>"
  )
set_target_properties(cpp_example PROPERTIES LINKER_LANGUAGE CXX)
add_test(cpp_example ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cpp_example)

# Compile Fortran example
if(ENABLE_FC_SUPPORT)
  add_executable(fortran_example
    fortran_example.f90
    ${PROJECT_SOURCE_DIR}/api/xcfun.F90
    )
  target_link_libraries(fortran_example
    XCFun::xcfun
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
