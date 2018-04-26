# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE
# Simplified and adapted by Roberto Di Remigio

#.rst:
#
# Enables code coverage and defines the list of appropriate flags for the GNU
# compiler.
#
# autocmake.yml configuration::
#
#   docopt: "--coverage Enable code coverage [default: OFF]."
#   define: "'-DENABLE_CODE_COVERAGE={0}'.format(arguments['--coverage'])"

option_with_print(ENABLE_CODE_COVERAGE "Enable code coverage" OFF)

if(ENABLE_CODE_COVERAGE)
  if(NOT CMAKE_BUILD_TYPE MATCHES "[Dd]ebug")
    message(WARNING "Code coverage analysis results with an optimized (non-Debug) build may be misleading")
  endif()

  if(NOT CMAKE_C_COMPILER_ID MATCHES GNU)
    message(FATAL_ERROR "Code coverage analysis only allowed with the GNU compiler collection!")
  endif()

  find_program(GCOV_PATH gcov)
  if(NOT GCOV_PATH)
    message(FATAL_ERROR "Code coverage analysis requires gcov!")
  endif()

  list(APPEND CODE_COVERAGE_FLAGS
    "-fprofile-arcs"
    "-ftest-coverage"
    )
endif()
