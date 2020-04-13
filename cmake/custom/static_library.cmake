#.rst:
#
# Enables creation of static library.
#
# Variables modified (provided the corresponding language is enabled)::
#
#   XCFun_CXX_FLAGS
#
# autocmake.yml configuration::
#
#   docopt: "--static Build as static library [default: False]."
#   define: "'-DBUILD_SHARED_LIBS={0}'.format(not arguments['--static'])"

option_with_print(BUILD_SHARED_LIBS "Build as shared library" ON)
