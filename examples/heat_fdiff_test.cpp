/*! \file heat_fdiff_test.cpp 
\brief A C++ program for testing the heat_fdiff library functions.

The program may be executed as follows to run all or individual tests.
 */
 
/*! \page heat_fdiff_test Test: Functions for testing heat_fdiff.cpp operation
\dontinclude heat_fdiff_test.cpp

\section intro Introduction
The program contained in heat_fdiff_test.cpp offers the user a series of tests for
the functionality of the critial heat_fdiff.cpp functions. The following section details
methods for implementing each or all of the tests.

\section build Building heat_fdiff_test
These testing procedures use the same makefile included in the examples directory. The 
function may be built with the following command: <tt>make heat_fdiff_test</tt>.

\section tests Test Implementation
The following list describes each of the available tests and a brief description. Note, that each test may be run sequentially by specifing a single processor or run with any number of processors, 3 is used here as an example.

-# <tt>mpirun -n 3 ./heat_fdiff_test all</tt> Implements all of the tests below.
-# <tt>mpirun -n 3 ./heat_fdiff_test get_coeff</tt> Tests that the b-vector coefficients are computed correctly
-# <tt>mpirun -n 3 ./heat_fdiff_test T</tt> Tests that the T-vector is being produced correctly with ghost values (only if the number of processes is 1 or 3 does this test produce data for comparision).
-# <tt>mpirun -n 3 ./heat_fdiff_test b</tt> Tests that the b-vector is updated correctly using T vector.
-# <tt>mpirun -n 3 ./heat_fdiff_test bvec_applyflux</tt> Tests that the flux is properly added to the b-vector.
*/ 
 
// Add standard files
#include <iostream.h>
#include <math.h> 
#include <time.h>
#include <string.h>

// Add the PETsc related files
#include "petscvec.h"	// PETsc vectors
#include "petscmat.h" 	// PETsc matrices
#include "petscksp.h" 	// PETsc KSP solver

// Add the solution and snow related files 
#include "heat_fdiff.h"	// Tools for solving 1-D heat equation
#include "SnowVar.h"	// Classes defining snow related variables
#include "Flux.h"		// Class for computing the various flux terms for a snowpack

// Define the help
static char help[] = "Test program for the heat_fdiff functions.\n";

// Begin the main program
int main(int argc, char *argv[])
{

// Initialize PETSc and define the process (r) and total number of processes (np)
int r, np;
PetscInitialize(&argc,&argv,(char *)0,help);
MPI_Comm_rank(MPI_COMM_WORLD, &r);
MPI_Comm_size(MPI_COMM_WORLD, &np);	
	
// get_coeff test
if (strcmpi(argv[1],"all") == 0 || strcmpi(argv[1],"get_coeff") == 0){
	
	// Indicate that test is starting
	printf("%s\n","TESTING: get_coeff");

	// Define necessary variables
	double c[4];
	
	// Run the function
	get_coeff(c, 100, 0.1, 1000, 1, 1);
	
	// Check the results
	if (c[0] == 0.1 &&  c[1] == 100000 && c[2] == 100000.1 && c[3] == 99999.9){
		printf("  %s\n\n","TEST PASSED: get_coeff");
		}
	else{
		printf("  %s\n\n","TEST FAILED: get_coeff");
	}
}

// createWithGhosts and assembleGhostVec (use "T" or "createWithGhosts" at command line)
if (strcmpi(argv[1],"all") == 0 || strcmpi(argv[1],"T") == 0 || strcmpi(argv[1],"createWithGhosts") == 0){
	
	// Indicate that test is starting
	if (r == 0){
		printf("%s\n","TESTING: createWithGhosts and assembleGhostVec");
	}
	
	// Cancel if number of processors is not 1 or 3
	if (r == 0 && np != 1 && np != 3){
		printf("!WARNING! Use 1 or 3 process to compare results with truth vectors !WANRING!\n");
	}
	
	// Initialize variables
	int i, low, high, n;	// Integer indices (i,low,high), local size of T (n)
	int N = np*4;			// Gloval size of T vector
	Vec T, TlocVec;			// storage location for global vec T and local vector TlocVec;
	double *Tloc;			// pointer for storing local values of T vector
	
	// Create the T vector: T[i] = -i
	createWithGhosts(&T,N);						// Calls heat_fdiff.cpp function
	VecGetLocalSize(T,&n);						// Size of local vector
	VecGetOwnershipRange(T, &low, &high);		// Global limits of the local vector
	
	// Loops through the vector and assigns values
	for (i = low; i<high; i++){					
		double val;	
		val = -i;
		VecSetValues(T,1,&i,&val,INSERT_VALUES);
	}
	
	// Assembles the vector
	assembleGhostVec(&T);						

	// Extract the local vector and store the values in Tloc
	VecGhostGetLocalForm(T,&TlocVec);			// Extracts the local values with ghosts to TlocVec
	VecGetArray(TlocVec,&Tloc);					// Extracts the local vector into an array, Tloc
	
	// Display the desired vector results
	if (r == 0 && np == 1){
		printf(" T vector should be: T = {0,-1,-2,-3, 0,-3}\n");
		printf(" The computed values are:\n");
	}
	else if (r == 0 && np == 3) {
		printf(" T vector should be for process [0]: T = {0,-1,-2,-3, 0,-4}\n");
		printf(" T vector should be for process [1]: T = {-4,-5,-6,-7, -3,-8}\n");
		printf(" T vector should be for process [2]: T = {-8,-9,-10,-11, -7,-11}\n");
		printf(" The computed values are:\n");
	}
	
	// Print the local forms of the T vector 
	printf("  Process [%d]\n",r);
	for (i=0; i<n+2; i++){	
		if (i < n){
			PetscSynchronizedPrintf(PETSC_COMM_SELF,"   T[%d] = %g\n",i,PetscRealPart(Tloc[i]));
		}
		else {
			PetscSynchronizedPrintf(PETSC_COMM_SELF,"   T[%d] = %g (ghost)\n",i,PetscRealPart(Tloc[i]));
		}
	}
	printf("\n");

	// Return the local array to the global temperature vector
	VecRestoreArray(TlocVec,&Tloc);			// completes use of local array Tloc
	VecGhostRestoreLocalForm(T,&TlocVec);	// completes the useage of loval vector TlocVec
	VecDestroy(&TlocVec);					// destroys the TlocVec
}

// bvec_update test
if (strcmpi(argv[1],"all") == 0 || strcmpi(argv[1],"bvec_update") == 0 || strcmpi(argv[1],"b") == 0 ){

	// Indicate that test is starting
	if (r == 0){
		printf("%s\n","TESTING: bvec_update");
	}

	// Define necessary variables
	int N = np*4;
	Vec T, bvec;
	
	// Create the necessary PETsc vectors
	createWithGhosts(&T,N);	// creates temperature vector (see heat_fdiff.h)	
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DETERMINE,N,&bvec);	// creates b vector of knowns
	
	// Set the initial temperature for all of the layers to a constant value	
	VecSet(T,-10);			// PETsc function for vector initialization
	assembleGhostVec(&T);	// function for assembling the ghosted vector (see heat_fdiff.h)

	// Run the function
	bvec_update(&bvec,&T, 1, 1, 1, -15);
	
	// Compare the results
	if (r == 0){
		printf("  %s\n","The desired results are: b = {-18,-20,-20,...,-20,-15}");
		printf("  %s\n","The actual results are:");
	}
	
	VecView(bvec,PETSC_VIEWER_STDOUT_WORLD);

	if (r == 0){
		printf("  %s\n"," ");
	}	
}
	
// bvec_applyflux test
if (strcmpi(argv[1],"all") == 0 || strcmpi(argv[1],"bvec_applyflux") == 0){
 	
	// Indicate that test is starting
	printf("%s\n","TESTING: bvec_applyflux");
 
	// Initizlize variables
	int N = np*4;	// Size of vector
	Vec u,v;		// PETSc vectors
 
	// Create and assemble a vector of all ones
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DETERMINE,N,&v);
	VecSet(v,1);
	VecAssemblyBegin(v);
	VecAssemblyEnd(v);
	
	// Copy this vector to u
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DETERMINE,N,&u);
	VecCopy(v,u);
	
	// Test the function (adds u and v together)
	bvec_applyflux(&v,&u);

	// Output results
	if (r == 0){
		printf("  The resulting vector should be a vector containing only 2's\n");
		printf("  The computed vector is:\n\n");
	}
	VecView(v,PETSC_VIEWER_STDOUT_WORLD);
}	

// Finialize the use of PETsc and MPI
	PetscFinalize();	
	return(0);
}




