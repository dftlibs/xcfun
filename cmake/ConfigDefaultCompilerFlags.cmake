if (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -fno-rtti -fno-exceptions ")
	set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_CXX_FLAGS_RELEASE "-g -O2 -DNDEBUG -Wno-unused")

	if (XCFUN_NO_STDC++)
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-threadsafe-statics")
	endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
	set (CMAKE_CXX_FLAGS "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177 -fno-rtti -fno-exceptions")
	set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
	set (CMAKE_CXX_FLAGS_RELEASE "-g -O2 -DNDEBUG")
endif ()

if (CMAKE_COMPILER_IS_GNUC)
	set (CMAKE_C_FLAGS "-Wall")
	set (CMAKE_C_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
	set (CMAKE_C_FLAGS_RELEASE "-g -O2 -DNDEBUG -Wno-unused")
elseif (CMAKE_C_COMPILER_ID MATCHES Intel)
	set (CMAKE_C_FLAGS "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
	set (CMAKE_C_FLAGS_DEBUG "-g -O0")
	set (CMAKE_C_FLAGS_RELEASE "-g -O2 -DNDEBUG")
endif ()

if (CMAKE_COMPILER_IS_GNUFORTRAN)
	set (CMAKE_Fortran_FLAGS "-Wall -Jfortran")
	set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
	set (CMAKE_Fortran_FLAGS_RELEASE "-g -O2")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
endif ()

# Take care of updating the cache for fresh configurations
if (NOT DEFINED DEFAULT_COMPILER_FLAGS_SET)
	mark_as_advanced(DEFAULT_COMPILER_FLAGS_SET)
	set_property(CACHE CMAKE_CXX_FLAGS 
		PROPERTY VALUE ${CMAKE_CXX_FLAGS})
	set_property(CACHE CMAKE_CXX_FLAGS_RELEASE 
		PROPERTY VALUE ${CMAKE_CXX_FLAGS_RELEASE})
	set_property(CACHE CMAKE_CXX_FLAGS_DEBUG 
		PROPERTY VALUE ${CMAKE_CXX_FLAGS_DEBUG})

	set_property(CACHE CMAKE_Fortran_FLAGS 
		PROPERTY VALUE ${CMAKE_Fortran_FLAGS})
	set_property(CACHE CMAKE_Fortran_FLAGS_RELEASE 
		PROPERTY VALUE ${CMAKE_Fortran_FLAGS_RELEASE})
	set_property(CACHE CMAKE_Fortran_FLAGS_DEBUG 
		PROPERTY VALUE ${CMAKE_Fortran_FLAGS_DEBUG})

	set (DEFAULT_COMPILER_FLAGS_SET ON 
		CACHE INTERNAL "Flag that the default compiler flags have been set.")
endif()
