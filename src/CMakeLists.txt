
include_directories(taylor functionals)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/include)

set (BASE_LIB_SRCS 
	xcfun.cpp
	evaluators.cpp
	functional.cpp
	parameters.cpp
	potential.cpp
	fortran.c
	)

file (GLOB FUNCTIONAL_SRCS functionals/*.cpp) 
set (GENERATED_FILES 
	${CMAKE_CURRENT_BINARY_DIR}/include/xcfun_autogen.h
	${CMAKE_CURRENT_BINARY_DIR}/Fortran/xcfun_autogen.f90
	)

add_library(xcfun ${BASE_LIB_SRCS}
	${FUNCTIONAL_SRCS} 
	${GENERATED_FILES})

add_executable(geninterface geninterface.cpp)

add_custom_target(interfaces ALL 
	DEPENDS 
	${CMAKE_CURRENT_BINARY_DIR}/include/xcfun_autogen.h
	${CMAKE_CURRENT_BINARY_DIR}/Fortran/xcfun_autogen.f90
)

add_custom_command (
	OUTPUT 
		${CMAKE_CURRENT_BINARY_DIR}/include/xcfun_autogen.h
		${CMAKE_CURRENT_BINARY_DIR}/Fortran/xcfun_autogen.f90
	DEPENDS 
		${CMAKE_CURRENT_BINARY_DIR}/include
		${CMAKE_CURRENT_BINARY_DIR}/Fortran
		geninterface 
		functionals/list_of_parameters.h 
		${PROJECT_SOURCE_DIR}/include/xcfun.h
	COMMAND geninterface
	)

add_custom_command (
	OUTPUT 
		${CMAKE_CURRENT_BINARY_DIR}/include
	COMMAND ${CMAKE_COMMAND} -E make_directory 
		${CMAKE_CURRENT_BINARY_DIR}/include
	)
add_custom_command (
    OUTPUT 
        ${CMAKE_CURRENT_BINARY_DIR}/Fortran
    COMMAND ${CMAKE_COMMAND} -E make_directory 
        ${CMAKE_CURRENT_BINARY_DIR}/Fortran
    )

install(TARGETS xcfun ARCHIVE DESTINATION lib)

