#.rst:
#
# autocmake.yml configuration::
#
#   docopt: "--enable-examples Enable code examples [default: False]."
#   define: "'-DENABLE_EXAMPLES={0}'.format(arguments['--enable-examples'])"

option(ENABLE_EXAMPLES "Enable code examples" OFF)

if(ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
