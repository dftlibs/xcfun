# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

#.rst:
#
# Enables code coverage by appending corresponding compiler flags.
#
# autocmake.yml configuration::
#
#   docopt: "--coverage Enable code coverage [default: OFF]."
#   define: "'-DENABLE_CODE_COVERAGE={0}'.format(arguments['--coverage'])"

option_with_print(
  NAME
    ENABLE_CODE_COVERAGE
  MESSAGE
    "Enable code coverage"
  DEFAULT
    OFF
  )

if(ENABLE_CODE_COVERAGE)
  if(NOT CMAKE_BUILD_TYPE STREQUAL "debug")
    message(WARNING "Code coverage analysis results with an optimized (non-Debug) build may be misleading")
  endif()

  find_program(GCOV_PATH gcov)
  if(NOT GCOV_PATH)
    message(FATAL_ERROR "Code coverage analysis requires gcov!")
  endif()
endif()
