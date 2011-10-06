/*! \file heat_fdiff.cpp 
\brief Source code for the functions designed to solve the 1-D heat equation with the finite difference method, details are provided in heat_fdiff.h.

Note: In general, these functions are passed pointers to the various vectors and matrices on which operations are being performed. This is done to minimize the potential for large vectors or matrices being copied and destoryed. As such, the * modifier is utilized. If comparing to the PETsc manual, this may be useful.
 */

// Add standard files
#include <iostream.h>
#include <math.h> 

// Add the PETsc related files
#include "petscvec.h"	// PETsc vectors
#include "petscmat.h" 	// PETsc matrices
#include "petscksp.h" 	// PETsc KSP solver

// Add the solution and snow related files
#include "heat_fdiff.h"	// Tools for solving 1-D heat equation

// Computes the a,b,c,d coefficents defined by Slaughter (2010)
int get_coeff(double *c, double rho, double k, double cp, double dz, double dt)
{
		c[0] = k/pow(dz,2);		// a
		c[1] = (rho * cp)/dt;	// b
		c[2] = c[1] + c[0];		// c
		c[3] = c[1] - c[0];		// d
		return(0);
}		

// Computes in parallel the vector of knowns using the ghosted temperature vector
int bvec_update(Vec *bvec, Vec *T, double a, double d, double qs, double Tbottom)
{
	int i, r, np;			// loop index, current process, total number of processes
	int low, high; 			// indices storing the range of global indices for a local vector
	int nlocal;				// size of local array
	int row = 0;			// index for inserting into PETsc vectors
	double T123[3];			// storage vector for temperatures {T(z-1), T(z), T(z+1)}
	double val;				// storage value for inserting into PETsc vectors
	double *Tloc;			// Pointer 
	
	// Define the current process and the total number of processes		
	MPI_Comm_rank(MPI_COMM_WORLD, &r);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	// Extract the local values of the T vector into the Tloc array
	VecGetLocalSize(*T,&nlocal);			// Size of local vector
	VecGetOwnershipRange(*T, &low, &high);	// Global limits of the local vector
	VecGetArray(*T,&Tloc);					// Extracts the local vector

	// Loop through each value of the local vector and compute the b vector component
	for (i=0; i<nlocal; i++) {
		row = low + i;	// Global index of the local component
	
		// Case for the first value of b, b[0]
		if(r == 0 && i == 0){						
			T123[1] = Tloc[0];
			T123[2] = Tloc[1];
			val = d*T123[1] + a*T123[2] + 2*qs; 	// b[0]
		} 
		
		// Case for the last value of b vector
		else if (r == (np-1) && i == (nlocal-1)){ 	
			val = Tbottom;							// b[N-1], N = global vector size
		}
		
		// General cases middle values of gloval b vector (b[1] to b[N-2])
		else { 										
			if (i == 0){				// case for first local entry
				T123[0] = Tloc[nlocal];		// 1st ghost (T[z-1])
				T123[1] = Tloc[i];			// T(z)
				T123[2] = Tloc[i+1];		// T[z+1]
			}
			else if (i == (nlocal-1)){	// case for middle local entries
				T123[0] = Tloc[i-1];		// T[z-1]
				T123[1] = Tloc[i];			// T[z]
				T123[2] = Tloc[nlocal+1];	// T[z+1]
			}
			else {						// case for last local entry
				T123[0] = Tloc[i-1];		// T[z-1]
				T123[1] = Tloc[i];			// T[z]
				T123[2] = Tloc[i+1];		// T[z+1]
			}
			val = a/2*T123[0] + d*T123[1] + a/2*T123[2]; // b value for general case
		}	
		
		// Insert the compute b vector value (all cases above)
		VecSetValues(*bvec,1,&row,&val,INSERT_VALUES);
	}
	
	// Return the local array to the global temperature vector
	VecRestoreArray(*T,&Tloc);
	
	// Assemble the b vector
	VecAssemblyBegin(*bvec);
	VecAssemblyEnd(*bvec);
	return(0);
}

// Adds the vector of absorbed heat flux to b vector of knowns
int bvec_applyflux(Vec *bvec, Vec *qabs)
{
	// Perform addition of vectors
	int a = 1;					// multiplier
	VecAXPY(*bvec, a, *qabs);	// y = y + a*x
	
	// Assemble the new b vector
	VecAssemblyBegin(*bvec);
	VecAssemblyEnd(*bvec);
	return(0);
}

// Assemble the global stiffness matrix in parallel
int Kmat_update(Mat *A,double a, double c)
{

	//  Define the variables
	int i, j;			// row loop index, column index,
	int r, N;		 	// current process, dimension of matrix (N x N_
	int low, high; 		// indices storing the range of global indices for a local vector
	double val;			// storage value for inserting into PETsc matrix object
	
	// Initize necessary variables
	MatGetSize(*A,&N,&N);			 		// global size of the stiffnexx matrix
	MPI_Comm_rank(MPI_COMM_WORLD, &r);		// the current process
	MatZeroEntries(*A);						// zero out the stiffness matrix
	MatGetOwnershipRange(*A, &low, &high);	// global index range of the local matrix

	// Loop through each row of the local stiffness matrix
	for(i = low; i < high; ++i) {
		// Case for the first row: A[0,0] and A[0,1]
		if(i == 0){
			MatSetValues(*A,1,&i,1,&i,&c,ADD_VALUES);	// A[0,0]
			j = 1;
			val = -a;
			MatSetValues(*A,1,&i,1,&j,&val,ADD_VALUES); // A[0,1]
		}
		
		// Case for the end of the matrix: A[N-1,N-1]
		else if(i == (N - 1)){
			val = 1.0;
			MatSetValues(*A,1,&i,1,&i,&val,ADD_VALUES);
	
		}
		
		// General case, diagonal terms
		else{
			MatSetValues(*A,1,&i,1,&i,&c,ADD_VALUES);	// A[i,i]
			j = i - 1;
			val = -a/2;
			MatSetValues(*A,1,&i,1,&j,&val,ADD_VALUES);	// A[i,i-1]
			j = i + 1;
			MatSetValues(*A,1,&i,1,&j,&val,ADD_VALUES); // A[i,i+1]
		}
	}
	
	// Assemble the stiffness matrix
	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd  (*A,MAT_FINAL_ASSEMBLY);
	return(0);
}

// Solves for the current temperature vector (x) given K and b: Ax = b 
int solve_temp(Vec *T, Vec *bvec, Mat *A)
{
	// Define variables
	int N;		// N = size of global vector
	Vec Tnew;	// PETsc vector for storing the computed temperatures (x)
	
	// Determine the global vector size
	VecGetSize(*T,&N);
	
	// Create the parallel Tnew vector
	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DETERMINE,N,&Tnew);
	
	// Define and solve the linear equations (see PETsc manual)
	KSP ksp; 							// Identifies KSP solution variable
	KSPCreate(PETSC_COMM_WORLD,&ksp);	// Creates the KSP solution object
	
	// Set tolerances
	KSPSetTolerances(ksp,PETSC_DEFAULT, PETSC_DEFAULT, 
		PETSC_DEFAULT, PETSC_DEFAULT); 
	KSPSetFromOptions(ksp);				
	
	// Assigns A matrix to solver
	KSPSetOperators(ksp,*A,*A,DIFFERENT_NONZERO_PATTERN); 
	
	// Solves the equation and outputs to Tnew
	KSPSolve(ksp,*bvec,Tnew); 				
	
	// Assign Tnew = T for the next time step
	VecCopy(Tnew,*T);
	
	// Destroy Tnew vector and zero out b vector and stiffness matrix
	VecDestroy(Tnew);
	VecZeroEntries(*bvec);
	MatZeroEntries(*A);
	return(0);
}

// Create a temperature vector with ghosted values
int createWithGhosts(Vec *x, int N)
{
	// Define variables
	int nlocal; 	// local size of vector
	int gidx[2]; 	// global indices for ghosts
	int r; 			// current process number
	int np; 		// total number of processes 
		
	// Determine process and total number of processors
	MPI_Comm_rank(MPI_COMM_WORLD, &r);	
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	// Return and error if the number of processes exceeds the length of vector
	if(np > N){
		SETERRQ(1,"The number of processors must not exceed the vectors size!");
	}

	// Calculate the desired number of local values per processes
	nlocal = ceil(double(N)/double(np));

	// Define the global indices for the ghosts
	// Case for first process
	if(r == 0){ 
		gidx[0] = 0; gidx[1] = nlocal;
	} 
	
	// Case for the last process
	else if(r == (np - 1)){
		gidx[0] = (r * nlocal) - 1; 
		gidx[1] = N - 1;
		nlocal = N - r * nlocal;
	}
	
	// Case for the middle processes
	else{
		gidx[0] = (r * nlocal) - 1;	
		gidx[1] = ((r + 1) * nlocal);	
	}	
	
	// Create the ghosted vector
	VecCreateGhost(PETSC_COMM_WORLD,nlocal,PETSC_DECIDE,2,gidx,x);
	return(0);
};

// Assembles the ghosted vector
int assembleVec(Vec *x)
{
	// Assembles the vector
	VecAssemblyBegin(*x);
	VecAssemblyEnd(*x);
	
	// Updates the ghost values
	VecGhostUpdateBegin(*x,INSERT_VALUES,SCATTER_FORWARD);
	VecGhostUpdateEnd(*x,INSERT_VALUES,SCATTER_FORWARD);
	
	return(0);
}

// Appends the temperature vector to a file
int output_temp(Vec *T, int init, char *filename, int n, double dz, int nt, double dt)
{
	// Define variables
	int r;			// current process
	double Tout;	// temperature to output
	int k = 0;		// loop index
	
	// Determine the current process
	MPI_Comm_rank(MPI_COMM_WORLD, &r);
	
	// Utilze the first process for creating the output file
	if( r == 0 ){
		// Define file and open the file for appending
		FILE * fid;
		fid = fopen(filename,"a");
		
		// Write the numerial solution information if the init flag is zero
		if(init == 0){
			fprintf(fid, "n = %i\n",n); 			// The number of layers 
			fprintf(fid, "dz = %f\n",dz); 			// The layer thickness (m)
			fprintf(fid, "nt = %i\n",nt);			// Number of time steps
			fprintf(fid, "dt = %f\n",dt);			// Time step (sec)
		}

		// Write the temperature data
		for(k=0; k<=n-1; k++){
			VecGetValues(*T,1,&k,&Tout);
			fprintf(fid,"%g\n",Tout);
		}
		fclose(fid);
	}
	return(0);
}