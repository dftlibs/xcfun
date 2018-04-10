option_with_default(
  NAME
    XCFun_XC_MAX_ORDER
  MESSAGE
    "Maximum order of derivatives of the exchange-correlation kernel"
  DEFAULT
    3
  )
# Make sure user selected a valuer larger than 2
if(DEFINED XCFun_XC_MAX_ORDER AND XCFun_XC_MAX_ORDER LESS 3)
  message(STATUS "${XCFun_XC_MAX_ORDER} not a valid value for maximum order of XC kernel derivatives! Resetting to its default value 3")
  set(XCFun_XC_MAX_ORDER 3 CACHE STRING "Maximum order of derivatives of the exchange-correlation kernel" FORCE)
endif()

# Hardcode to share, rather than use CMAKE_INSTALL_DATAROOTDIR as the latter
# might resolve to a place not recognized by CMake
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")

if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
endif()
file(TO_NATIVE_PATH "${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}/xcfun" PYMOD_INSTALL_FULLDIR)
