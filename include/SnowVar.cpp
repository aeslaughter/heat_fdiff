/*! \file SnowVar.cpp 
\brief Source code for the classes used to define the constant values needed for the finite difference solution for snow as defined and documented in the header file SnowVar.h.
 */
 
#include "SnowVar.h" 			// Associated header file
//! \todo Make these one class with subclasses; add member functions to allow variables to change with time
// Initilize the default values for the snow, environmental, and numerical constants.
SnowConstants::SnowConstants(){		
	Ls = 2833; 		
	Ke = 0.0023; 	
	Kh = Ke; 		
	MvMa = 0.622; 	
	Ra = 0.287; 		
	Rv = 0.462; 		
	T0 = -5; 		
	e0 = 0.402; 
	k  = 0.1;
	p = 150;
	cp = 2000;
	emis = 0.98;
	};
	
EnvirConditions::EnvirConditions(){	
	qsw = 500;	
	qlw = 200;	
	alpha = 0.9;
	kappa = 100;	
	Ta = -10;
	Patm = 102;	
	cpa = 1003;
	Vw = 5;			
	Tbot = -10;		
	Tstart = -10;
	RH = 50;
	};
	
FiniteDiffParam::FiniteDiffParam(){	
	n  = 10; 			
	dz = 0.01; 		
	dt = 60;			
	nt = 60;			
	};
	
	


