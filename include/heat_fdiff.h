/*! 
\mainpage Parallel 1-D Finite Difference Solution Library with Application for Layered Snowpacks.

\section intro_sec Introduction
 The programs included here are some basic functions for solving the heat equation in 1-D using the finite difference method.  Additionall, included are functions for applying the method to model a layered snowpack, as presented by  <A href="http://proquest.umi.com/pqdlink?Ver=1&Exp=10-01-2016&FMT=7&DID=2048139271&RQT=309&attempt=1&cfc=1"> Slaughter (2010)</A>. The code presented here offers a simplified solution, that is the input does not vary with time. However, this code could easily be modified to make some or all variables dependant on time.
 
 The main purpose for developing these programs was primarily to learn C++, parrallel programming, and documentation using 
 <A HREF="http://www.stack.nl/~dimitri/doxygen/index.html">Doxygen</A>. As such, the code is likely far more complicated and commented than need be for such a simple problem.  Additionally, the author hopes that these codes will be useful for others learning C++ and parallel programming at the same time.
 
 \section files Associated Source Files
 The functions contained in heat_fdiff.cpp were designed to be somewhat generic in nature and could be applied to solve the 1-D heat equation in general.  That being said, the functions were setup based on the derivation presented by  <A href="http://proquest.umi.com/pqdlink?Ver=1&Exp=10-01-2016&FMT=7&DID=2048139271&RQT=309&attempt=1&cfc=1">Slaughter (2010)</A>, which was a specific solution for snow.
 
 The classes and associated member functions contained in Flux.cpp and SnowVar.cpp are specific to the application of the heat equation to a layered snowpack.
 
 The main.cpp file contains an example of how to use this code with specific application to a layered snowpack.

 \section eqn Governing Equation
The governing differential equation which is solved numerically using the finite difference is as follows:
 
 \f$ \rho c_p \frac{\partial T}{\partial t} = k  \frac{\partial^2 T}{\partial z^2} -  \frac{\partial q}{\partial z}\f$,
 
 where  \f$ \rho \f$ is density,  \f$ c_p \f$ is specific heat,  \f$ k \f$ is thermal conductivity,  \f$ T \f$ is temperature,  \f$ t \f$ is time, and  \f$ z \f$ is depth. Complete details of the numerical solution is given in <A href="http://proquest.umi.com/pqdlink?Ver=1&Exp=10-01-2016&FMT=7&DID=2048139271&RQT=309&attempt=1&cfc=1">Chapter 5 of Slaughter (2010)</A>.
 
 \section Installation
 To compile the library the makefile in the main directory should be executed. The user should change the path to this library in the makefile, this is indicated in the makefile comments. The paths PETsc and MPI should not require any changes if the library is being run on the MPDC cluster.

In the examples folder contains an example that demonstrates how to use this library to solve for the temperatures of a layered snowpack. This code is compiled using the makefile in the examples directory. Again, the path to this library should be editted.  Details explaining the example are given here: \ref example1. 
 
 \section tools Additional Tools
 Included in the main directory of the library is a MATLAB file, plotT.m, that reads and plots the outputed data files from the finite difference library.
 
 The documentation was created with the Doxygen Wizard GUI. The settings file associated is included in the main directory:  doc_settings.
 
 \section References
Slaughter, A.E. <A href="http://proquest.umi.com/pqdlink?Ver=1&Exp=10-01-2016&FMT=7&DID=2048139271&RQT=309&attempt=1&cfc=1">"Numerical analysis of conditions necessary for near-surface snow metamorphism."</A>
Ph.D. Dissertation, Montana State University, 2010, 604 pages.
 */
 
/*! \file heat_fdiff.h
\brief Header file for heat_fdiff.cpp with code for solving the 1-D heat equatin in parallel.

The functions contained here are intented to be generic in nature. That is, in the example provided these functions are used to solve for the temperatures in a layered snowpack, but as written these functions could be applied to any number of heat related problems.
 */
 
 /*! \fn get_coeff(double *, double, double, double, double, double);
\brief Adds absorbed flux vector to the b vector of knowns.

<B>SYNTAX:</B><ul> 
	<li>get_coeff(c, rho, k, cp, dz, dt);</ul>

<B>INPUT:</B><ul>
	<li> c = Pointer to a an array with four elements
	<li> rho = density
	<li> k = thermal conductivity
	<li> cp = specific heat capacity
	<li> dz = Layer thickness
	<li> dt = Time step</ul>
	
These coefficients are referred to by Slaughter (2010) as a, b, c, and d. Here they are contained in order in a single array c.	
*/
int get_coeff(double *, double, double, double, double, double);


/*! \fn bvec_update(Vec *, Vec *, double, double, double, double);
\brief Computes the vector of knowns without the absorbed heat flux

<B>SYNTAX:</B><ul> 
	<li>bvec_update(&b,&T, c[0], c[3], qs, Tstart);</ul>

<B>INPUT:</B><ul>
	<li> &b = Pointer to a PETsc vector object that contains the knowns (see bvec_update)
	<li> c[0] = The first, a in Slaughter (2010), coefficient (see get_coeff)
	<li> c[3] = The fourth, d in Slaughter (2010), coefficient (see get_coeff)
	<li> qs = The heat flux at the surface node
	<li> Tstart = The initial staring temperature (assumed constant for the entire system)</ul>
*/
int bvec_update(Vec *, Vec *, double, double, double, double);


/*! 
\fn bvec_applyflux(Vec *, Vec *);
\brief Adds absorbed flux vector to the b vector of knowns.

<B>SYNTAX:</B><ul> 
	<li>bvec_applyflux(&b, &q);</ul>

<B>INPUT:</B><ul>
	<li> &b = Pointer to a PETsc vector object that contains the knowns (see bvec_update)
	<li> &q = Pointer to a PETsc vector object containing absorbed heat to be added</ul>
*/
int bvec_applyflux(Vec *, Vec *);


/*! \fn Kmat_update(Mat *,double, double);
\brief Builds or updates the stiffness matrix

<B>SYNTAX:</B><ul> 
	<li>Kmat_update(&A,a,b);</ul>

<B>INPUT:</B><ul>
	<li> &A = Pointer to a PETsc matrix object containing the stiffness matrix (see Kmat_update)
	<li> a,b = Coefficients used in stiffness and b-matrix, see get_coeff and Slaughter (2010) for detals</ul>	
*/
int Kmat_update(Mat *,double, double);


/*! \fn solve_temp(Vec *, Vec *, Mat *);
\brief Solves the linear equation, Ax = b, for x in parallel.

<B>SYNTAX:</B><ul> 
	<li>solve_temp(&T, &b, &A);</ul>

<B>INPUT:</B><ul>
	<li> &T = Pointer to a PETsc vector object that will store the computed temperatures
	<li> &b = Pointer to a PETsc vector object that contains the knowns (see bvec_update)
	<li> &A = Pointer to a PETsc matrix object containing the stiffness matrix (see Kmat_update)</ul>	
*/
int solve_temp(Vec *, Vec *, Mat *);


/*! \fn createWithGhosts(Vec *, int);
\brief Creates a parallel vector with padding or ghosted values.

<B>SYNTAX:</B><ul> 
	<li>createWithGhosts(&T), where T is a PETsc vector object.</ul>

<B>INPUT:</B><ul>
	<li> &T = Pointer to a PETsc vector object</ul>
		
<B>DESCRIPTION:</B>\n
The finite difference solution for the transient 1-D heat requires that the temperature vector for the current time step be used to compute the future temperature vector. With the exception at the boundaries, this calculation states that the future temperature at time t + 1 at node i, T(i,t+1), is a function of the of the three other temperatures at the current time: T(i-1,t), T(i,t), and T(i+1,t). As such a parallel distributed vectors must have overlapping values for proper computation.


For example, if the vector [0,1,2,3,4,5,6,7,8] was parsed to three processors, then using the following would be produced with this function:

1: [0,1,2|,0,3]\n
2: [3,4,5|,2,6]\n
3: [6,7,8|,5,8]\n

The vertical line is used to distiguish between the local vector and the ghost values for this vector, which are the values before and after the vector. At the boundaries the first and last entries are simply repeated, but not used in calculations. Note, the function bvec_update utilizes these padded vectors, this function only provides instructions for assembly and the function assembleVec actually performs the assembly.
*/
int createWithGhosts(Vec *, int);


/*! \fn assembleVec(Vec *);
\brief Assembles a ghosted PETsc vector created with the function createWithGhosts

<B>SYNTAX:</B><ul> 
	<li>assembleVec(&T), where T is a PETsc vector object.</ul>

<B>INPUT:</B><ul>
	<li> &T = Pointer to a PETsc vector object</ul>	
*/
int assembleVec(Vec *);


/*! \fn output_temp(Vec *, int, char *, int, double,int, double);
\brief Outputs a temperature vector to a file.

<B>SYNTAX:</B><ul> 
	<li>output_temp(&T, trig, filename, dz, n, dt, nt);</ul>

<B>INPUT:</B><ul>
	<li> &T = Pointer to a PETsc vector object
	<li> trig = 0 or 1, if zero the header dz, n, dt, and nt are printed
	<li> dz = The distance between nodes
	<li> n = The number of nodes
	<li> dt = Time step in seconds
	<li> nt = Number of times steps</ul>	
*/
int output_temp(Vec *, int, char *, int, double, int, double);
