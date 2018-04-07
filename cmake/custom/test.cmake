option_with_print(
  NAME
    ENABLE_TESTALL
  MESSAGE
    "Enable compilation of testall"
  DEFAULT
    ON
  )

if(ENABLE_TESTALL)
  include(CTest)
  enable_testing()
  add_subdirectory(test)
endif()
