/*! \file main.cpp 
\brief A 1-D Finite Difference Solution for Snow.

 The following includes a 1-D finite difference solution of the heat equation including short-wave radiation penetration. The following is an example of how to utilize the classes and functions defined in the library.
 */
 
// MAIN PROGRAM PREAMBLE
// Add standard files
#include <iostream.h>
#include <math.h> 
#include <time.h>

// Add the PETsc related files
#include "petscvec.h"	// PETsc vectors
#include "petscmat.h" 	// PETsc matrices
#include "petscksp.h" 	// PETsc KSP solver

// Add the solution and snow related files
#include "heat_fdiff.h"	// Tools for solving 1-D heat equation
#include "SnowVar.h"			// Classes defining snow related variables
#include "Flux.h"				// Class for computing the various flux terms for a snowpack

// Define the help
static char help[] = "Solves, in parallel, the 1-D heat equation for a layed snowpack using finite difference.\n";
//using namespace std; // Uncomment this to use cout, etc. directly, else std::cout is required.

// BEGIN THE MAIN PROGRAM
int main(int argc, char *argv[]){

// Define the start time
	time_t starttime, endtime;
	time(&starttime);
	
// Define the output filename. 
	/* The default is output.dat, but this may be changed via the command line. The first command line input, if it exists, is assumed to be the desired filename. For example:
	[user]$ ./main results.txt.
	*/
	char *filename;				// Pointer to a character array
	if(argv[1] != NULL){ 		// Case without a command line input
		filename = argv[1];
	}
	else{						// Case with a command line filename input
		filename = "output.dat";
	}

// Create objects for apply solution to a layered snowpack
	FiniteDiffParam D;	// class that stores parameters related to the numerical solution
	EnvirConditions E;	// class containing variables related to enviornmental conditions
	SnowConstants C;	// class containing variables related to snow material properties
	Flux F;				// class with member functions for computing heat flux terms
	
// Define the variables to be used
	int i, t;			// indices used for loops
	int low, high; 		// indices storing the range of global indices for a local vector
	int N = (D.n) + 1; 	// the total number of layers/elements used
	int row;			// generic index for passing information to PETsc objects
	double val;			// generic variable used for passing values to PETsc objects
	double c[4];		// array for storing the a,b,c,d coefficients (see Slaughter 2010)
	double qs, Ts;		// surface flux and surface temperature
	Vec T, bvec, qabs;	// PETsc vectors for temperature, vector of knowns, and absorbed flux
	Mat K;				// PETsc matrix for the stiffness matrix
		
// Initializes PETsc, including MPI
	PetscInitialize(&argc,&argv,(char *)0,help);	

// Create the necessary PETscvectors
	createWithGhosts(&T,N);	// creates temperature vector (see finite_diff_tools.h)
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DETERMINE,N,&bvec);	// creates b vector of knowns
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DETERMINE,N,&qabs); // creates vector of absorbed flux
	
// Create the stiffnes matrix (see PETsc manual for additional details) 
	MatCreateMPIAIJ(PETSC_COMM_WORLD, PETSC_DETERMINE, PETSC_DETERMINE, N, N, 
		3, PETSC_NULL, 0, PETSC_NULL, &K); 
		
// Set the initial temperature for all of the layers to a constant value	
	VecSet(T,E.Tstart);	// PETsc function for vector initialization
	assembleVec(&T);	// function for assembling the ghosted vector (see finite_diff_tools.h)
	
// Output the initial temperature to the desired file
	output_temp(&T, 0, filename, D.n, D.dz, D.nt, D.dt);
	
// Begin computation of new temperatures by stepping through time
	for(t=0; t <= D.nt; t++)
	{
	// Construct the absorbed heat flux vector
		/* Note, both the absorbed (qabs) and surface (qs) variables are in units of W/m^3. The true heat flux values are divided by the layer thickness in the Flux member functions. Thus, although incorrect, the term heat flux is used here for convience.*/
		VecGetOwnershipRange(qabs, &low, &high);	// global limits of the local vector
		for(i = low; i < high; ++i) {				// loop over the gather limits
			val = F.qabsorbed(D,E,C,i);				// compute the absorbed flux at this layer
			VecSetValues(qabs, 1, &i, &val, INSERT_VALUES);	// Insert the value into the vector
		}
		
	// Assemble the absrobed heat flux vector
		VecAssemblyBegin(qabs);
		VecAssemblyEnd(qabs);
		
	// Comppute the surface heat flux
		row = 0;					// Index used to extract a value from a PETsc vector
		VecGetValues(T,1,&row,&Ts); // Extracts the surface temperature, Ts
		qs = F.qsurface(D,E,C,Ts);	// Calculate the surface heat flux
		
	// Compute the a,b,c,d coefficients (see Slaughter 2010)
		/* Note, as is the case for the stiffness matrix and portions of the b vector, these values are not a function of time, thus could be outside of the loop. They are placed here for future modification to make this example and library more robust.*/
		get_coeff(c, C.p, C.k, C.cp, D.dz, D.dt); 

	// Updates b vector with the coefficients, surface flux, and fixed boundary temperature
		bvec_update(&bvec,&T, c[0], c[3], qs, E.Tstart);
		
	// Adds the abaorbed heat flux to the b vector
		bvec_applyflux(&bvec,&qabs);
		
	// Creates the stiffness matrix
		Kmat_update(&K,c[0],c[2]);
	
	// Solves the Ax = b equation and returns the new temperature to T vector
		solve_temp(&T, &bvec, &K);

	// Writes the new temperature vector to the file
		output_temp(&T, 1, filename, D.n, D.dz, D.nt, D.dt);		
	};
	
// Finialize the use of PETsc and MPI
	PetscFinalize();				

// Define the end time
	time(&endtime);
	printf("Execution time: %f.\n", difftime(endtime,starttime));    
	return(0);
	
}





