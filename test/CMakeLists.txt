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
    ${PROJECT_SOURCE_DIR}/api/xcfun.f90
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

if(XCFun_ENABLE_PYTHON_INTERFACE)
  list(APPEND _pytest_files
    ${CMAKE_CURRENT_SOURCE_DIR}/test_xcfun.py
    )
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/tests)
  foreach(_file IN LISTS _pytest_files)
    get_filename_component(__file ${_file} NAME)
    add_custom_target(symlink_${__file} ALL
      COMMAND
        ${CMAKE_COMMAND} -E create_symlink
        ${_file}
        ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/tests/${__file}
      )
  endforeach()

  add_test(
    NAME python-interface
    COMMAND ${PYTHON_EXECUTABLE} -m pytest -x -rws -vv ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/tests
    WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/..
    )

  install(
    FILES
      "${_pytest_files}"
    DESTINATION
      ${PYMOD_INSTALL_FULLDIR}/tests
    PERMISSIONS
      OWNER_EXECUTE
      OWNER_WRITE
      OWNER_READ
      GROUP_EXECUTE
      GROUP_READ
    )
endif()

