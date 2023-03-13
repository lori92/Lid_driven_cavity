# -*- mode: makefile -*-

#  This sample Makefile can be used to compile PETSc applications
#  Copy this file to your source directory as "Makefile" and modify as needed.
#  See also $PETSC_DIR/share/petsc/Makefile.user for the preferred approach
#  You must set the environmental variable(s) PETSC_DIR (and PETSC_ARCH if PETSc was not configured with the --prefix option)
##
MAX_THREADS=4

PETSC_DIR =/root/petsc
PATH_LIBS=/root/libs
PATH_OPENBLASLIB=/root/libs_openBlas/lib
PATH_OPENBLASLINC=/root/libs_openBlas/include

#   You can set specific values below but that should rarely need to
CFLAGS		 =O3
FFLAGS		 =
CPPFLAGS         =
FPPFLAGS         =

CXX        = g++
CXXFLAGS   = -Wall -fbounds-check -std=c++14
#CXXFLAGS   = -MMD -MP -I$(PATH_OPENBLASLINC) -pthread -fopenmp -O3 -funroll-all-loops -fexpensive-optimizations -ftree-vectorize -fprefetch-loop-arrays -floop-parallelize-all -ftree-parallelize-loops=$(MAX_THREADS) -m64 -c -Wall#
STRIP      = strip
LDFLAGS    = 


EXECUTABLE = cavity

app : LDLIBS += -lstdc++ 



#$(EXECUTABLE): $(OBJS)
#	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)
#	$(STRIP) $(EXECUTABLE)



SRCS=main.cpp 
OBJS=main.o mesh.o 

#$(EXECUTABLE): $(OBJS)
#	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDLIBS) -L$(PATH_OPENBLASLIBS) -lopenblas -lgfortran

#$(EXECUTABLE): $(OBJS)
#	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDLIBS) -L$(PATH_LIBS) -L$(PATH_LIBS) -llapack -lblas -lgfortran

#$(OBJS): $(SRCS)
#	$(CXX) $(CPPFLAGS) -c  $(SRCS) $(LDLIBS) 

#mesh.o: mesh.cpp
#	$(CXX) $(CPPFLAGS) -c  $(SRCS) $(LDLIBS) 

#	g++ -c mesh.cpp                # translates frac.cpp into object code, frac.o 
#	g++ -c main.cpp                # translates main.cpp into object code, main.o



#mesh.o: mesh.cpp
#	$(CXX) $(CPPFLAGS) -c  mesh.cpp $(LDLIBS) 


#main.o: main.cpp
#	$(CXX) $(CPPFLAGS) -c  main.cpp $(LDLIBS)	

#$(EXECUTABLE): $(OBJS)
#	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDLIBS) -L$(PATH_LIBS) -L$(PATH_LIBS) -llapack -lblas -lgfortran

SOURCES = main.cpp write_output.cpp mesh.cpp  predictor_step_new.cpp poisson_solver.cpp finiteVolume.cpp numerics.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: clean all
.DEFAULT_GOAL := all

all: cavity

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ -c

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LDLIBS) -L$(PATH_LIBS) -L$(PATH_LIBS) -llapack -lblas -lgfortran

clean:
	rm -f *.o $(EXECUTABLE) *.d
    
-include $(OBJECTS:.o=.d)


#  When linking in a multi-files with Fortran source files a.F90, b.c, and c.cxx
#  You may need to use
#
# app : a.o b.o c.o
# 	$(LINK.F) -o $@ $^ $(LDLIBS)

# If the file c.cxx needs to link with a C++ standard library -lstdc++ , then
# you'll need to add it explicitly.  It can go in the rule above or be added to
# a target-specific variable by uncommenting the line below.
#include ${PETSC_DIR}/lib/petsc/conf/variables

#  To access the PETSc variables for the build, including compilers, compiler flags, libraries etc but
#  manage the build rules yourself (rarely needed) comment out the next lines
#include ${PETSC_DIR}/lib/petsc/conf/rules
#include ${PETSC_DIR}/lib/petsc/conf/test



