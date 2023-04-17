# -*- mode: makefile -*-

#  This sample Makefile can be used to compile PETSc applications
#  Copy this file to your source directory as "Makefile" and modify as needed.
#  See also $PETSC_DIR/share/petsc/Makefile.user for the preferred approach
#  You must set the environmental variable(s) PETSC_DIR (and PETSC_ARCH if PETSc was not configured with the --prefix option)
##
BLAS_LAPL_LIBSS=/root/libs
PETSC_DIR=/root/petsc

# Include Petsc-defined variables

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test


#   You can set specific values below but that should rarely need to
CFLAGS		 =O3
FFLAGS		 =
CPPFLAGS         =
FPPFLAGS         =


CXX_STD = -std=c++11
CXX_CFLAGS = ${CXX_STD} ${CXX_FLAGS} ${PETSC_CCPPFLAGS} 
LIBS = ${PETSC_LIB} ${BLAS_LAPL_LIBS}
CXX_LFLAGS = ${CXX_STD}

STRIP      = strip
LDFLAGS    = 


#EXECUTABLE = cavity

#app : LDLIBS += -lstdc++ 

SOURCES =   FVmesh.cpp FVIO.cpp FVmatrix.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: default allclean



all: prog
prog: ${OBJECTS} FVmain.o 
	@echo "--- COMPILE FVmain  with ${CXX}  ---"
	$(CXX) ${OBJECTS} FVmain.o -o prog ${LIBS} ${CXX_CFLAGS}

FVmain.o: ${OBJECTS} FVmain.cpp
	@echo "--- COMPILE FVmain  with ${OBJECTS}  ---"
	${CXX} -o FVmain.o -c FVmain.cpp ${LIBS} ${CXX_CFLAGS}
	@echo "==============="

FVmatrix.o: FVmesh.o FVmatrix.cpp
	@echo "--- COMPILE FVmatrix---"
	${CXX} -o FVmatrix.o -c FVmatrix.cpp ${CXX_CFLAGS}
	@echo "==============="

FVIO.o: FVmesh.o FVIO.cpp
	@echo "--- COMPILE FVIO ---"
	${CXX} -o  FVIO.o -c FVIO.cpp ${CXX_CFLAGS}
	@echo "==============="

FVmesh.o: FVmesh.cpp
	@echo "--- COMPILE MESH---"
	${CXX} -o FVmesh.o -c FVmesh.cpp ${CXX_CFLAGS}
	@echo "==============="


#numerics.o: mesh.o numerics.cpp
#	@echo "--- COMPILE NUMERICS---"
#	${CXX} -o  numerics.o -c numerics.cpp ${LIBS} ${CXX_CFLAGS}
#	@echo "==============="

#predictor_step_new.o: mesh.o  predictor_step_new.cpp
#	@echo "--- COMPILE NUMERICS---"
#	${CXX} -o  predictor_step_new.o -c predictor_step_new.cpp ${CXX_CFLAGS}
#	@echo "==============="

#poisson_solver.o: mesh.o  poisson_solver.cpp
#	@echo "--- COMPILE poisson_solver---"
#	${CXX} -o  poisson_solver.o -c poisson_solver.cpp ${CXX_CFLAGS}
#	@echo "==============="



#finiteVolume.o: mesh.o finiteVolume.cpp
#	@echo "--- COMPILE finiteVolume ---"
#	${CXX} -o  finiteVolume.o -c finiteVolume.cpp ${CXX_CFLAGS}
#	@echo "==============="

#BoundaryConditions.o: mesh.o BoundaryConditions.cpp
#	@echo "--- COMPILE BoundaryConditions ---"
#	${CXX} -o  BoundaryConditions.o -c BoundaryConditions.cpp ${CXX_CFLAGS}
#	@echo "==============="




