add_executable(testall
  testall.cpp
  )
target_link_libraries(testall
  xcfun
  )
target_compile_options(testall
  PRIVATE
    "${XCFun_CXX_FLAGS}"
    "$<$<CONFIG:Debug>:${XCFun_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${XCFun_CXX_FLAGS_RELEASE}>"
    "$<$<BOOL:${ENABLE_CODE_COVERAGE}>:${XCFun_CXX_FLAGS_COVERAGE}>"
  )
add_test(
  NAME testall
  COMMAND $<TARGET_FILE:testall>
  )

if(XCFUN_PYTHON_INTERFACE)
  list(APPEND _pytest_files
    test_xcfun.py
    )
  foreach(_file IN LISTS _pytest_files)
    add_custom_target(symlink_${_file} ALL
      COMMAND
        ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/${_file}
        ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../${_file}
      )
  endforeach()

  add_test(
    NAME python-interface
    COMMAND ${PYTHON_EXECUTABLE} -m pytest -x -rws -vv $<TARGET_FILE_DIR:_xcfun>/..
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

