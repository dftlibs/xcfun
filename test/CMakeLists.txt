
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/src
    )
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(funeval funeval.c)
target_link_libraries(funeval xcfun m)

add_executable(testall testall.c)
target_link_libraries(testall xcfun m)

if (XCFUN_ENABLE_TESTS)
	add_test(xcfun testall)
endif()
