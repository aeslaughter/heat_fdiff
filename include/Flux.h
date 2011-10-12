/*! \class Flux
	\brief A class for gathering and computing the various flux terms for associated with the snowpack.
	
The purpose of the Flux class is to provided the necessary tools to (1) compute the surface flux for application of the Nuemann boundary condition (constant flux) and (2) to compute the absorbed short-wave radiation flux for a given layer of the snowpack. This calculations are performed in the member functions qsurface and qabsorbed, respectively. All other member functions are used by these two functions.
*/

class Flux{
	public:
	double pa;  //!< Air density (kg/m^3); see qsurface
	
	/*! A  class member that uses the Clausius-Clayperon equation for computing the water-vapor pressure over a surface.
	
	<B>SYNTAX:</B><ul> 
		<li>ei = clausius(C,Ti)</ul>
	
	<B>INPUT:</B><ul> 
		<li> C = SnowConstants object
		<li>Ti (double) = the temperture (in degrees C) at which the vapor-pressure shall be computed.</ul>
	
	<B>INPUT:</B><ul> 
		<li>ei (double) = the vapor-pressure (in kPa) computed.</ul>
	*/
	double clausius(SnowConstants, double);
	
	/*! A class member function for setting the value for the ensity of air
	
	<B>SYNTAX:</B><ul> 
		<li>airdensity(E,C)</ul>
	
	<B>INPUT:</B><ul>
		<li> E = EnvirConditions object
		<li> C = SnowConstants object</ul>
		
	<B>OUTPUT:</B><ul>
		<li> pa = The density of air (kg/m^3)</ul>
	
	*/
	void airdensity(EnvirConditions, SnowConstants);
	
	/*! A class member function for the latent heat flux
	
	<B>SYNTAX:</B><ul> 
		<li>qe = latent(E,C,Ts)</ul>
	
	<B>INPUT:</B><ul>
		<li> E = EnvirConditions object
		<li> C = SnowConstants object
		<li> Ts = Snow surface temperature (C)</ul>
		
	<B>OUTPUT:</B><ul>
		<li> qe = Latent heat flux (W/m^2)</ul>
	*/
	double latent(EnvirConditions, SnowConstants, double);
	
	/*! A class member function for the sensible heat flux
	
	<B>SYNTAX:</B><ul>
		<li>qh = sensible(E,C,Ts)</ul>
	
	<B>INPUT:</B><ul>
		<li> E = EnvirConditions object
		<li> C = SnowConstants object
		<li> Ts = Snow surface temperature (C)</ul>
		
	<B>OUTPUT:</B><ul>
		<li> qh = Sensible heat flux (W/m^2)</ul>
	*/
	double sensible(EnvirConditions, SnowConstants, double);
	
	/*! A class member function for the longwave heat flux
	
	<B>SYNTAX:</B><ul>
		<li>qlw = longwave(E,C,Ts), where F is an instance of the Flux class.</ul>
	
	<B>INPUT:</B><ul>
		<li> E = EnvirConditions object
		<li> C = SnowConstants object
		<li> Ts = Snow surface temperature (C)</ul>
		
	<B>OUTPUT:</B><ul>
		<li> qlw = Longwave heat flux (W/m^2)</ul>
	*/
	double longwave(EnvirConditions, SnowConstants, double);
	
	/*! A class member for building the absorbed short-wave flux vector.
		
	There is no ouput for this member function, but the function uses the *q pointer to update the absorbed flux vector.
		
	<B>SYNTAX:</B><ul> 
		<li>qa = F.qabsorbed(D,E,C,i), where F is an instance of the Flux class.</ul>
	
	<B>INPUT:</B><ul>
		<li> D = FiniteDiffParam object
		<li> E = EnvirConditions
		<li> C = SnowConstants object
		<li> i = Layer number
		</ul>
		
	<B>OUTPUT:</B><ul>
		<li> qa = The absorbed flux for layer i (W/m^3)</ul>
	*/	
		
	double qabsorbed(FiniteDiffParam, EnvirConditions, SnowConstants, int);
		
	/*! A class member for building the surface flux vector.
		
	<B>SYNTAX:</B><ul> 
		<li>qs = F.qsurface(D,E,C,Ts), where F is an instance of the Flux class.</ul>
	
	<B>INPUT:</B><ul>
		<li> D = FiniteDiffParam object
		<li> E = EnvirConditions
		<li> C = SnowConstants object
		<li> Ts = Snow surface temperature (C)</ul>
			
	<B>OUTPUT:</B><ul>
		<li> qs = Surface heat flux (W/m^3)</ul>
	
	*/
	double qsurface(FiniteDiffParam, EnvirConditions, SnowConstants, double);
};

 