# List options
option_with_print(
  NAME
    ENABLE_64BIT_INTEGERS
  MESSAGE
    "Enable 64-bit integers"
  DEFAULT
    OFF
  )

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
# FIXME Localize use of this option
add_definitions("-DXC_MAX_ORDER=${XCFun_XC_MAX_ORDER}")
