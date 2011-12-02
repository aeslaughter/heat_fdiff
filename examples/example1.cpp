/*! \page example1 Example 1: Layered Snowpack
\dontinclude example1.cpp

\section ex1_intro Introduction
The following includes a 1-D finite difference solution of the heat equation including short-wave radiation penetration. The following is an example of how to utilize the classes and functions defined in the library. The complete source code file is provided here: example1.cpp.  

\section ex1_desc Description of Code
\subsection ex1_preamble Preamble
First the necessary header files are loaded, including the necessary standard libraries. Then the necessary PETsc related headers are included as perscribed in the PETsc documentations for each of the functions used. Finally, the headers associated with this library are included, where heat_fdiff.h provides the general functions and the others are used for direct application to snow.

\snippet example1.cpp preamble

\subsection ex1_main Begin Main Function and User Input
The main function begins with initizining a time variable for documenting the computation time. Then the command line inputs are handled. By default the temperature data produced is saved in output.dat, but the use may include the desired name from the command line as an additional argument (e.g., ./main data.txt).
\snippet example1.cpp main

\subsection ex1_var Define Objects and Variables
Four objects are created, the first three, D, E, and C, are basically storage structures for the necessary parameters need to solve the heat equation for snow (see SnowVar.h). The last object, F, is a Flux object that contains member functions for computing the various flux terms for the snow problem (see Flux.h).

Next, all of the variables used throughout the program are defined, including the necessary PETsc vectors and matrix for performing parallel calculations.
\snippet example1.cpp var

\subsection ex1_build Initilize the PETsc Objects
The PETsc vectors---T, qabs, and bvec---and the stiffness matrix K, must be initilzed. This is accomplished be first initilizing the PETsc package itself, which also initializes MPI. 

The temperature vector has special requirements, thus the function from this library, createWithGhosts, was used (see heat_fdiff.h). The other vectors and K matrix are created in standard fashion with PETsc functions. The stiffness matrix is diagnol sparse matrix, well suited for the MatCreateMPIAIJ function provided with PETsc.

Finally, the temperature vector is set to the constant starting temperature and assembled using the assebleVec function from this library (see heat_fdiff.h).
\snippet example1.cpp build

\subsection ex_output1 Initial Output
Before begining the time loop, the initial temperature vectors is output using the output function from this library (see heat_fdiff.h).
\snippet example1.cpp out1

\subsection ex_loop Begin Time Loop
The temperature for each time step is solved within a loop.
\snippet example1.cpp loop

\subsection ex_snow Compute the Snow Related Terms
The two flux terms, the surface flux qs, and the absorbed flux, qabs are computed as they pertain to a layered snow pack using the Flux object (see Flux.h) supplied with this library. First the absorbed flux is computed for each layer, in parallel, with the qabsorbed member function of the Flux object. This vector is then assembled with PETsc's assembly commands. Finally the surface flux is computed with the qsurface member function of the Flux class. 
\snippet example1.cpp snow

\subsection ex_solve Solve for New Temperatures
Seven steps, in this example are required to solve and record the new temperatures. 
	-# THe vector and matrix coefficients are computed based on the material properties and the numerical time step and layer thickness.
	-# Based on these coefficients and the flux, the vector of knowns (bvec) is updated.
	-# The absorbed flux vector (qabs) is added to the vector of knowns (bvec).
	-# The stiffness matrix, K, is updated.
	-# The linear system is solved using the solve_temp function of this library (see heat_fdiff.h).
	-# The new temperature is outputed to the desired file.
	-# Finally, the temperature array is adjusted not to allow the snow to get above freezing.
	
\snippet example1.cpp solve

\subsection ex_finish Complete the Program
After completing all of the desired time steps within the loop the program is completed by finalize the PETsc package and outputing the computation time.
\snippet example1.cpp finish

\subsection ex_results Example Results
Running this example as is produces the output included in example1.dat. Then using the plotT.m MATLAB function the following contour plot is created from this file, which displays the temperature as a function of depth and time is created.
\image html example1.png

*/ 

/*! \example example1.cpp 
\brief An example solution for the 1-D finite difference libraries applied to snow. See \ref example1 for details.
*/

/*! \file example1.cpp 
\brief An example solution for the 1-D finite difference libraries applied to snow. See \ref example1 for details.
*/

//! [preamble]
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
#include "SnowVar.h"	// Classes defining snow related variables
#include "Flux.h"		// Class for computing the various flux terms for a snowpack

// Define the help
static char help[] = "Solves, in parallel, the 1-D heat equation for a layed snowpack using finite difference.\n";
//using namespace std; // Uncomment this to use cout, etc. directly, else std::cout is required.
//! [preamble]

//! [main]
// BEGIN THE MAIN PROGRAM
int main(int argc, char *argv[])
{

// Define the start time
	time_t starttime, endtime;
	time(&starttime);

// Define the output filename. 
	char *filename;				// Pointer to a character array
	filename = "output.dat";

//! [main]

//! [var]		
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
//! [var]	

//! [build]
// Initializes PETsc, including MPI
	PetscInitialize(&argc,&argv,(char *)0,help);	

// Create the necessary PETsc vectors
	createWithGhosts(T,N);	// creates temperature vector (see heat_fdiff.h)	
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DETERMINE,N,&bvec);	// creates b vector of knowns
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DETERMINE,N,&qabs); // creates vector of absorbed flux
	
// Create the stiffnes matrix (see PETsc manual for additional details) 
	MatCreateMPIAIJ(PETSC_COMM_WORLD, PETSC_DETERMINE, PETSC_DETERMINE, N, N, 
		3, PETSC_NULL, 0, PETSC_NULL, &K); 
		
// Set the initial temperature for all of the layers to a constant value	
	VecSet(T,E.Tstart);	// PETsc function for vector initialization
	assembleGhostVec(T);	// function for assembling the ghosted vector (see heat_fdiff.h)
//! [build]	
	
//! [out1]	
// Output the initial temperature to the desired file
	output_temp(T, 0, filename, D.n, D.dz, D.nt, D.dt);
//! [out1]	
	
//! [loop]	
D.nt = 0;
// Begin computation of new temperatures by stepping through time
	for(t=0; t <= D.nt; t++)
	{
//! [loop]

//! [snow]	
	// Construct the absorbed heat flux vector
		/* Note, both the absorbed (qabs) and surface (qs) variables are in units of W/m^3. The true heat flux values are divided by the layer thickness in the Flux member functions. Thus, although incorrect, the term heat flux is used here for convience.*/
		VecGetOwnershipRange(qabs, &low, &high);			// global limits of the local vector
		for(i = low; i < high; ++i) {						// loop over the gather limits
			val = F.qabsorbed(D,E,C,i);						// compute the absorbed flux at this layer
			VecSetValues(qabs, 1, &i, &val, INSERT_VALUES);	// insert the value into the vector
		}
		
	// Assemble the absrobed heat flux vector
		VecAssemblyBegin(qabs);
		VecAssemblyEnd(qabs);
		
	// Comppute the surface heat flux
		row = 0;					// Index used to extract a value from a PETsc vector
		VecGetValues(T,1,&row,&Ts); // Extracts the surface temperature, Ts
		qs = F.qsurface(D,E,C,Ts);	// Calculate the surface heat flux
//! [snow]	

//! [solve]
	// Compute the a,b,c,d coefficients (see Slaughter 2010)
		/* Note, as is the case for the stiffness matrix and portions of the b vector, these values are not a function of time, thus could be outside of the loop. They are placed here for future modification to make this example and library more robust.*/
		get_coeff(c, C.p, C.k, C.cp, D.dz, D.dt); 

	// Updates b vector with the coefficients, surface flux, and fixed boundary temperature
		bvec_update(bvec, T, c[0], c[3], qs, E.Tbot);
		
	// Adds the abaorbed heat flux to the b vector
		bvec_applyflux(bvec,qabs);
	
	// Creates the stiffness matrix
		Kmat_update(K,c[0],c[2]);
	
	// Solves the Ax = b equation and returns the new temperature to T vector
		solve_temp(T, bvec, K);
		
	// Do not allow the snow to melt
		VecGetOwnershipRange(T, &low, &high);	// global limits of the local vector
		for(i = low; i < high; ++i) {			// loop over the gather limits
			VecGetValues(T,1,&row,&val); 		// extracts the surface temperature, Ts
			if(val > 0){
				val = 0;
				VecSetValues(T, 1, &i, &val, INSERT_VALUES);	// insert the value into the vector
			}
		}	
		assembleGhostVec(T); // assemble T vector
	
	// Writes the new temperature vector to the file
		output_temp(T, 1, filename, D.n, D.dz, D.nt, D.dt);	
//! [solve]

//! [finish]		
	}; // Complete the time loop
	
// Finialize the use of PETsc and MPI
	PetscFinalize();				

// Define the end time
	time(&endtime);
	printf("Execution time: %f.\n", difftime(endtime, starttime));    
	return(0);
}
//! [finish]