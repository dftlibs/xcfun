#This file probably requires GNU Make, but should be easily
#convertible to regular Make format.

# Compile time options:
# -DNDEBUG Turn off run time checks in production calculations to gain speed
# -DNO_STDCXX Do not use any part of the C++ standard 
#             library, to allow linking without it.
# -DWITH_QD Enable quad double precision functions. This requires
#           a _patched_ version of the QD library
# -DWITH_SINGLE Enable single precision interface. This doubles the
#               size of the binary, so be careful.

# C++ compiler and flags
CXX ?= g++
CXXFLAGS ?= -g -Wall -O #-DWITH_QD
CXXFLAGS+=-Iinclude -Isrc/taylor -Isrc/functionals -Llib -fno-rtti -fno-exceptions
LIBS= #-lqd

# F90 compiler and flags, used for example code
FC=gfortran
FFLAGS=-Wall -Llib -Jfortran # -DWITH_QD

xcfun_test: src/xc_fun.cpp $(wildcard src/functionals/*.h) include/xc_fun.h src/taylor/taylor.h
	$(CXX) $(CXXFLAGS)  src/xc_fun.cpp -o $@ $(LIBS)

lib: lib/libxc_fun.a

lib/libxc_fun.a: src/xc_fun.o $(wildcard src/functionals/*.h) include/xc_fun.h src/taylor/taylor.h
	$(CXX) $(CXXFLAGS) -c src/xc_fun.cpp -o src/xc_fun_lib.o $(LIBS) -DXCFUN_LIB	
	ar r $@ src/xc_fun_lib.o


stripped: lib/libxc_fun.a
	strip --strip-unneeded $<

funeval: test/funeval.cpp lib/libxc_fun.a
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS) -lxc_fun

fortran/example.o: fortran/xc_fun_module.o

fortran_example: fortran/example.o fortran/xc_fun_module.o lib/libxc_fun.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) -lstdc++

clean:
	rm -f `find . -name '*.o' -o -name '*~' -o -name '*.a'` fortran_example benchmark testxc


.SUFFIXES: .F90

.F90.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.F90
