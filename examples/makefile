# Define the compiler and compiler options
CC = g++
CFLAGS	= -g

# Define the paths for using PETsc and MPI 
include ${PETSC_DIR}/conf/variables

# Define path of the heat equation finite difference library (location of this file)
FDPATH = /cygdrive/c/users/slaughter/Documents/Cornell/software/heat_fdiff

# Define the path for the BOOST header files
BOOSTPATH = /cygdrive/c/cygwin/packages/boost_1_47_0

# Define the include directories for  this library
FDINCPATH = $(FDPATH)/include

# Define the list of paths and flags to pass to the compiler
INCPATH	= -I$(FDINCPATH) -I${BOOSTPATH} ${PETSC_CC_INCLUDES}

# Define the list of library to pass to the compiler
LIBPATH = ${PETSC_LIB}

# Define path to object and lib folders for use in compiler commands
OBJ = $(FDPATH)/object
LIB = $(FDPATH)/lib

# Define the paths for make to automatically search
VPATH = $(FDINCPATH) $(OBJ) $(LIB)

example1: example1.o
	$(CC) -o $@ example1.o $(OBJ)/*.o $(LIBPATH)	
	
example1.o: example1.cpp libheatfdiff.so
	$(CC) -c $(CFLAGS) $(INCPATH) $(INCFLAG) -I$(OBJ) example1.cpp	
	
heat_fdiff_test: heat_fdiff_test.o
	$(CC) -o $@ heat_fdiff_test.o $(OBJ)/*.o $(LIBPATH)	
	
heat_fdiff_test.o: heat_fdiff_test.cpp libheatfdiff.so
	$(CC) -c $(CFLAGS) $(INCPATH) $(INCFLAG) -I$(OBJ) heat_fdiff_test.cpp		
	
clean:
	rm -f example1 example1.o example1.exe example1.dat heat_fdiff_test heat_fdiff_test.o heat_fdiff_test.exe
	