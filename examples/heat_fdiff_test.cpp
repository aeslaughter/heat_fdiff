/*! \file heat_fdiff_test.cpp 
\brief A C++ program for testing the heat_fdiff library functions.

The program may be executed as follows to run all or individual tests.
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

// BEGIN THE MAIN PROGRAM
int main(int argc, char *argv[])
{

// GET_COEFF TEST
if (strcmpi(argv[1],"all") == 0 || strcmpi(argv[1],"get_coeff") == 0){
	// Define necessary variables
	double c[4];
	
	// Indicate that test is starting
	printf("%s\n","TESTING: get_coeff");
	
	// Run the function
	get_coeff(c, 100, 0.1, 1000, 1, 1);
	
	// Check the results
	if (c[0] == 0.1 &&  c[1] == 100000 && c[2] == 100000.1 && c[3] == 99999.9){
		printf("  %s\n","TEST PASSED: get_coeff");
		}
	else{
		printf("  %s\n","TEST FAILED: get_coeff");
	}
}

// BVEC_UPDATE TEST
if (strcmpi(argv[1],"all") == 0 || strcmpi(argv[1],"bvec_update") == 0){
	// Define necessary variables
//	double c[4];
	
	// Indicate that test is starting
	printf("%s\n","TESTING: bvec_update");
	
	// Run the function
//	get_coeff(c, 100, 0.1, 1000, 1, 1);
	
	// Check the results
/*	if (c[0] == 0.1 &&  c[1] == 100000 && c[2] == 100000.1 && c[3] == 99999.9){
		printf("  %s\n","TEST PASSED: bvec_update");
		}
	else{
		printf("  %s\n","TEST FAILED: bvec_update");
	}
*/	
}






return(0);
}




