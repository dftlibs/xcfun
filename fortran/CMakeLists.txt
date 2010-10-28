
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(example_fortran 
	example.f90
	xcfun_module.f90
	${PROJECT_BINARY_DIR}/src/fortran/xcfun_autogen.f90
	)

target_link_libraries(example_fortran xcfun m)
