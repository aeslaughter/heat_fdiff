# Define the compiler and compiler options
CC = g++
CFLAGS	= -w -O2

# Define the paths for using PETsc and MPI 
PETSCPATH = /software/petsc/petsc-2.3.3
PETSCARCH = linux-gnu
MPICHPATH = /opt/mpich-ch_p4-gcc-1.2.7

# EDIT THIS PATH!
# Define path of the heat equation finite difference library (location of this file)
FDPATH = /home/abh76/programs/heat_fdiff

# Define the include directories for PETsc and this library
PETSCINCPATH =  $(PETSCPATH)/include
PETSCARCHINCPATH = $(PETSCPATH)/bmake/$(PETSCARCH)
MPICHINCPATH = $(MPICHPATH)/include
FDINCPATH = $(FDPATH)/include

# Define the list of paths and flags to pass to the compiler
INCPATH	= -I$(FDINCPATH) -I$(PETSCINCPATH) -I$(PETSCARCHINCPATH) -I$(MPICHINCPATH) -I$(METISINCPATH) 
INCFLAG = -DPETSC_USE_DEBUG -DPETSC_USE_LOG -DPETSC_USE_BOPT_O

# Define the library paths
PETSCLIBPATH = -L$(PETSCPATH)/lib/$(PETSCARCH) -lpetscksp -lpetscmat -lpetscvec -lpetsc -lpetscdm
MPICHLIBPATH = -L$(MPICHPATH)/lib -lmpich 
FDLIBPATH = -L$(FDPATH)/lib -lheatfdiff

# Define the list of library to pass to the compiler
LIBPATH = $(PETSCLIBPATH) $(MPICHLIBPATH)  

# Define path to object and lib folders for use in compiler commands
OBJ = $(FDPATH)/object
LIB = $(FDPATH)/lib

VPATH = $(FDINCPATH) $(OBJ) $(LIB)	

example1: example1.o
	$(CC) -O -o $@ example1.o $(OBJ)/*.o $(LIBPATH)	
	
example1.o: example1.cpp libheatfdiff.so
	$(CC) -c $(CFLAGS) $(INCPATH) $(INCFLAG) -I$(OBJ) example1.cpp		
	
clean:
	rm -f example1 example1.o example2 example2.o
	