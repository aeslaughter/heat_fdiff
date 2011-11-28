# Define the compiler and compiler options
CC = g++
CFLAGS	= -w -O2 -g

# Define the paths for using PETsc and MPI 
include ${PETSC_DIR}/conf/variables

# Define path of the heat equation finite difference library (location of this file)
ifeq ($(PETSC_ARCH), win32-gnu)
	FDPATH = /cygdrive/c/users/slaughter/Documents/Cornell/software/heat_fdiff
else
	FDPATH = /home/andrew/Documents/win-documents/Cornell/software/heat_fdiff
endif

# Define the include directories for  this library
FDINCPATH = $(FDPATH)/include

# Define the list of paths and flags to pass to the compiler
INCPATH	= -I$(FDINCPATH) ${PETSC_CC_INCLUDES}

# Define the list of library to pass to the compiler
LIBPATH = ${PETSC_LIB}

# Define path to object and lib folders for use in compiler commands
OBJ = $(FDPATH)/object
LIB = $(FDPATH)/lib

# Define the paths for make to automatically search
VPATH = $(FDINCPATH) $(OBJ) $(LIB)

# Compile commands for the library file and each of the associated objects
libheatfdiff.so: heat_fdiff.o Flux.o SnowVar.o
	$(CC) -shared -o $(LIB)/$@ $(OBJ)/*.o $(LIBPATH)	

Flux.o: Flux.cpp SnowVar.o
	$(CC) -fPIC -c  $(CFLAGS) $(INCPATH) $(INCFLAG) $(FDINCPATH)/Flux.cpp -o $(OBJ)/$@
	
SnowVar.o: SnowVar.cpp
	$(CC) -fPIC -c $(CFLAGS) $(INCPATH) $(INCFLAG) $(FDINCPATH)/SnowVar.cpp -o $(OBJ)/$@
	
heat_fdiff.o: heat_fdiff.cpp
	$(CC) -fPIC -c $(CFLAGS) $(INCPATH) $(INCFLAG) $(FDINCPATH)/heat_fdiff.cpp -o $(OBJ)/$@
	
clean:
	rm -f $(OBJ)/*.o $(LIB)/*.so
	
