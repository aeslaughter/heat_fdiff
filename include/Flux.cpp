/*! \file Flux.cpp 
\brief Source code for Flux class used for finite difference solution for snow as defined and documented in the header file Flux.h.
 */
 
#include <math.h>		// Mathmatical functions library
#include "SnowVar.h"	// Snow Variables header file (required)
#include "Flux.h"		// Associated header file

using namespace std;

double Flux::clausius(SnowConstants C, double Ti)
	{
		double ei;
		ei = C.e0 * exp(C.Ls/C.Rv * (1/(C.T0+273.15) - 1/(Ti+273.15)));
		return(ei);
	};
	
void Flux::airdensity(EnvirConditions E, SnowConstants C)
	{
		pa = E.Patm/(C.Ra*(E.Ta+273.15)); 
	};
	
double Flux::latent(EnvirConditions E, SnowConstants C, double Ts)
	{
		double qe;  //<! Latent heat flux (W/m^2)
		//NOTE: The air density must be initilized before latent may be used
		qe = 1000*(C.MvMa*pa*C.Ls*C.Ke*E.Vw)/E.Patm * ((clausius(C,E.Ta)*E.RH/100) - clausius(C,Ts));
		return(qe);
	};
		
double Flux::sensible(EnvirConditions E, SnowConstants C, double Ts)
	{
		double qh;  //<! Sensible heat flux (W/m^2)
		//NOTE: The air density must be initilized before latent may be used
		qh = pa*E.cpa*C.Kh*E.Vw*(E.Ta - Ts); // Sensible heat flux		
		return(qh);
	};
	
double Flux::longwave(EnvirConditions E, SnowConstants C, double Ts)
	{
		double qlw;  //<!Longwave heat flux (W/m^2)
		qlw = E.qlw - C.emis * 5.670e-8 * pow((Ts+273.15),4); // Long-wave radiation
		return(qlw);
	};
				
double Flux::qabsorbed(FiniteDiffParam F, EnvirConditions E, SnowConstants C, int i)
	{
		double q1, q2, qabs;	//<! Intermediate and outputed flux terms
		q1 = (E.qsw * (1 - E.alpha)) * exp(-E.kappa * (i) * F.dz);
		q2 = (E.qsw * (1 - E.alpha)) * exp(-E.kappa * (i+1) * F.dz);
		qabs = (q1 - q2)/(F.dz);
		return(qabs);
	}; 

double Flux::qsurface(FiniteDiffParam F, EnvirConditions E, SnowConstants C, double Ts)
	{
		// NOTE: The qs scalar is initilized in the Flux class definition
		double qs;
		airdensity(E,C); // Initilizes the air density

		qs = (longwave(E,C,Ts) + latent(E,C,Ts) + sensible(E,C,Ts))/F.dz; // Total surface heat flux
		return(qs);
	};
	
			