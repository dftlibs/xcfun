
link_directories(${PROJECT_BINARY_DIR}/src)

# Dummy target to satisfy dependency
add_custom_command (
	OUTPUT 
	${PROJECT_BINARY_DIR}/src/fortran/xcfun_autogen.f90
	DEPENDS 
	${PROJECT_BINARY_DIR}/src/geninterface 
	)

add_executable(example_fortran 
	example.f90
	xcfun_module.f90
	${PROJECT_BINARY_DIR}/src/fortran/xcfun_autogen.f90
	)

set_target_properties(example_fortran 
	PROPERTIES LINKER_LANGUAGE Fortran
	)

target_link_libraries(example_fortran xcfun m)
