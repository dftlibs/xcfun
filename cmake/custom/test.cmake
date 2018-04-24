option_with_print(ENABLE_TESTALL "Enable compilation of testall" ON)

if(ENABLE_TESTALL)
  include(CTest)
  enable_testing()
  add_subdirectory(test)
endif()
