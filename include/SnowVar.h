/*! \file SnowVar.h 
\brief Header file for SnowVar.cpp.

Application of the heat equation to a layered snowpack requires the application of imperical relationships to describe the surface and absorbed heat flux.  As such, a number of parameters are required. These three classes served to help organize these parameters. As implemented these class are simple storage structures that do not vary with time.  However, these classes could easly be modified to allow any or all of the parameters to vary with time. 
 */
 
 /*! \class SnowConstants 
	\brief Class containing constants for snow.

	SnowConstants is a class that contains the various constants that are required for using the finite difference code for a layered snowpack. In general, these values do not need to be editted. This class is basically a data structure containing the constant values used for modeling snow.
 */
class SnowConstants{
	public:
		double Ls;	 //!< Latent heat of sublimation (kJ/kg)
		double Ke; 	 //!< Transfer coefficient for water vapor
		double Kh;	 //!< Transfer coefficient
		double MvMa; //!< Ratio of dry-air to water-vapor molecular weights
		double Ra;	 //!< Gas constant for air (kJ/kg/K)
		double Rv;	 //!< Gas constant for water-vapor (kJ/kg/K)
		double T0;	 //!< Reference temperature for vapor pressure (C)
		double e0;	 //!< Reference vapor pressure (kPa)
		double k; 	 //!< Thermal conductivity (W/m/K)
		double p; 	 //!< Snow density (kg/m^3)
		double cp;	 //!< Specific heat capacity (kJ/kg/K)
		double emis; //!< Emissivity for snow
		
	/*! SnowConstants constructor.

	When an SnowConstants object is created the constructor initilizes the variables with the values listed.
	*/	
	SnowConstants();
};
	
 /*! \class EnvirConditions
	\brief Class containing environmental variables for modeling snow.

	EnvirConditions is a class that contains the various environmental parameters that remain constant with time and depth. These parameters are initlized with values, thus offer quick access for a demo.  But, usually these terms are changed to meet the specific conditions the user is desiring to model. Similar to SnowConstants, this class is simply acting as a data structure storing the environmental conditions.
*/
class EnvirConditions{
	public:
		double qsw;		//!< Incoming short-wave radiation (W/m^2)
		double qlw;		//!< Incoming long-wave radiation (W/m^2)
		double alpha; 	//!< Snow albedo
		double kappa;	//!< Snow extinction coefficient (1/m)
		double Ta;		//!< Air temperature
		double Patm;	//!< Air pressure (kPa)
		double cpa;		//!< Specific heat capacity of air (J/kg/K)
		double Vw;		//!< Wind speed (m/s)
		double Tbot;	//!< Lower boundary condition (C)
		double Tstart;	//!< Snow starting temperature (C)
		double RH;		//!< Relative humidity (%)

	/*! EnvirConditions constructor.

		When an EnvirConditions object is created the constructor initilizes the variables with the values listed.
	*/
	EnvirConditions();
};

/*! \class FiniteDiffParam
	\brief Class containing parameters pertitnent to the finite difference method as well as the output file name.

	FiniteDiffParam is a class that contains the parameters need to use the finite difference method to solve the 1-D heat eqution progroslm. These parameters are initlized with values, thus offer quick access for a demo.  But, usually these terms are changed to meet the specific grid refinement that the user desires to model. Similar to SnowConstants and EnvirConditions, this class is simply acting as a data structure.
*/ 
class FiniteDiffParam{
	public:
		int n; 			//!< The number of layers 
		double dz; 		//!< The layer thickness (m)
		double dt;		//!< Time step (sec)
		int nt;			//!< Number of time steps to compute
		char output[];  //!< Output filename
			
	/*! FiniteDiffParam constructor.

	When an FiniteDiffParm object is created the constructor initilizes the variables with the values listed.
	*/
	FiniteDiffParam();
 };
 
