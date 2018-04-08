# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE
# Simplified and adapted by Roberto Di Remigio

#.rst:
#
# Enables OpenMP support.
#
# Variables used::
#
#   ENABLE_OPENMP
#   OPENMP_FOUND
#
# autocmake.yml configuration::
#
#   docopt: "--omp Enable OpenMP parallelization [default: False]."
#   define: "'-DENABLE_OPENMP={0}'.format(arguments['--omp'])"

option_with_print(
  NAME
    ENABLE_OPENMP
  MESSAGE
    "Enable OpenMP parallelization"
  DEFAULT
    OFF
  )

if(ENABLE_OPENMP)
  if(NOT OPENMP_FOUND)
    find_package(OpenMP REQUIRED)
  endif()
endif()
