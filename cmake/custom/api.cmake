option_with_print(
  NAME
    ENABLE_Fortran_INTERFACE
  MESSAGE
    "Enable Fortran interface"
  DEFAULT
    ON
  )
if(ENABLE_Fortran_INTERFACE)
  message(STATUS "Fortran interface ENABLED")
  add_subdirectory(fortran)
endif()

option_with_print(
  NAME
    ENABLE_PYTHON_INTERFACE
  MESSAGE
    "Enable Python interface"
  DEFAULT
    OFF
  )
if(ENABLE_PYTHON_INTERFACE)
  message(STATUS "Python interface ENABLED")
  add_subdirectory(python)
endif()
