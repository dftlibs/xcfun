
include_directories(
    taylor 
    functionals
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_BINARY_DIR}/src
    )

set (BASE_LIB_SRCS 
	xcfun.cpp
	evaluators.cpp
	functional.cpp
	parameters.cpp
	potential.cpp
	fortran.c
	)

file (GLOB FUNCTIONAL_SRCS functionals/*.cpp) 

add_library(xcfun ${BASE_LIB_SRCS}
	${FUNCTIONAL_SRCS} 
	${GENERATED_FILES}
	xcfun_autogen.h
	xcfun_autogen.f90
    )

add_executable(geninterface geninterface.cpp)

add_custom_target(interfaces ALL 
	DEPENDS 
	${CMAKE_CURRENT_BINARY_DIR}/xcfun_autogen.h
	${CMAKE_CURRENT_BINARY_DIR}/xcfun_autogen.f90
)

add_custom_command (
	OUTPUT 
		${CMAKE_CURRENT_BINARY_DIR}/xcfun_autogen.h
		${CMAKE_CURRENT_BINARY_DIR}/xcfun_autogen.f90
	DEPENDS 
		geninterface 
		functionals/list_of_parameters.h 
		${PROJECT_SOURCE_DIR}/include/xcfun.h
	COMMAND geninterface
	)

install(TARGETS xcfun ARCHIVE DESTINATION lib)

