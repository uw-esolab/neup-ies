/*************************************** FIT Water (NREL v2) ***************************************

Temperature Range: 273.2 K to 1,500.0 K
Pressure Range: 1.0 Pa to 50.0 MPa
FIT Version: 4cd1a62f8548

----------------------------------------------------------------------------------------------------
Copyright (c) 2016, Northland Numerics LLC
All rights reserved.

Use of this software in source and binary forms, with or without modification, is permitted for
Alliance for Sustainable Energy, LLC (the “Licensee”) and for third parties that receive this
software, with or without modification, directly from Licensee.  Licensee is permitted to
redistribute this software in source and binary forms, with or without modification, provided
that redistributions of source code must retain the above copyright notice and reservation of
rights, this paragraph, and the following disclaimer of warranties and liability.  Absent
separate written agreement with Northland Numerics LLC, redistribution of this software is not
permitted for third parties receiving this software, with or without modification, from Licensee
or from any other entity or individual.

THIS SOFTWARE IS PROVIDED BY NORTHLAND NUMERICS LLC "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL ANY WARRANTY BE CREATED IN
CONNECTION WITH THE SALE OF THIS SOFTWARE, UNLESS THE WARRANTY WAS CREATED SOLELY AND EXPRESSLY
IN A SEPARATE WRITTEN AGREEMENT SIGNED BY NORTHLAND NUMERICS LLC.  IN NO EVENT SHALL ANY WARRANTY
BE IMPUTED OR PRESUMED.  IN NO EVENT SHALL NORTHLAND NUMERICS LLC OR CONTRIBUTORS TO THIS SOFTWARE
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, BUSINESS INTERRUPTION; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
OR LOSS OF USE, DATA, OR PROFITS), HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE), ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

----------------------------------------------------------------------------------------------------

The file water_properties.c contains a number of functions that return water properties
calculated using an interpolated Helmholtz free energy and its analytical derivatives.  The main
thermodynamic property functions are:

  int water_TD( double T, double D, water_state * state )
  int water_TP( double T, double P, water_state * state )
  int water_PH( double P, double H, water_state * state )
  int water_PS( double P, double S, water_state * state )
  int water_HS( double H, double S, water_state * state )
  int water_TQ( double T, double Q, water_state * state )
  int water_PQ( double P, double Q, water_state * state )

The first two arguments correspond to the known independent properties:

  T -- temperature (K)
  D -- density (kg/m3)
  P -- pressure (kPa)
  H -- enthalpy (kJ/kg)
  S -- entropy (kJ/kg-K)
  Q -- quality on a mass basis

The property functions require a pointer to a water_state structure, which is updated with the
thermodynamic properties at the state defined by the independent variables.  If an error occurs, 
the function returns a value other than 0.  More information about the specific return value is
provided by the function water_error_message.

Additional functions that are available (all inputs and return values are type double):

  water_visc(D, T) -- viscosity (uPa-s) as a function of density and temperature
  water_cond(D, T) -- thermal conductivity (W/m-K) as a function of density and temperature

Warning: The above functions will return -9.0e99 if the input is not valid.

Notes:
  1) The thermodynamic state is only explicitly defined in temperature and density.  Therefore, any
     other combination of known properties requires iteration and may result in a state that does
     not exactly correspond to the provided properties.  For this reason, the values that are
     set for the specified properties may not be identical to the inputs.  However, the values do
     exactly correspond to the thermodynamic state defined by the returned temperature and density.
  2) The get_water_info function provides various fluid-specific and FIT-specific information by
     setting the values of the water_info structure provided to it.

For more information, contact Northland Numerics at: support@nnumerics.com

***************************************************************************************************/

#ifndef FIT_WATER_H_INCLUDED
#define FIT_WATER_H_INCLUDED

typedef struct water_info
    {
    double molar_mass;        // molar mass of fluid (kg/kmol)
    double T_critical;        // critical temperature (K)
    double D_critical;        // critical density (kg/m3)
    double P_critical;        // critical pressure (kPa)
    double temp_lower_limit;  // lowest available temperature (K)
    double temp_upper_limit;  // highest available temperature (K)
    double pres_lower_limit;  // lowest available pressure (kPa)
    double pres_upper_limit;  // highest available pressure (kPa)
    double sat_temp_min;      // minimum saturation temperature (K)
    double sat_pres_min;      // minimum saturation pressure (kPa)
    }
    water_info;

typedef struct water_state
    {
    double temp;          // temperature (K)
    double pres;          // pressure (kPa)
    double dens;          // density (kg/m3)
    double qual;          // quality (-)
    double inte;          // internal energy (kJ/kg)
    double enth;          // enthalpy (kJ/kg)
    double entr;          // entropy (kJ/kg-K)
    double cv;            // specific heat at const. volume (kJ/kg-K)
    double cp;            // specific heat at const. pressure (kJ/kg-K)
    double ssnd;          // speed of sound in fluid (m/s)
    double sat_vap_dens;  // saturated vapor density (kg/m3)
    double sat_liq_dens;  // saturated liquid density (kg/m3)
    }
    water_state;

// Information functions.
const char * water_error_message( int error_code );
void get_water_info( water_info * info );

// Thermodynamic property functions. 
// T: Temperature; D: density; P: Pressure; H: enthalpy; S: entropy; Q: quality; 
int water_TD( double T, double D, water_state * state ); 
int water_TP( double T, double P, water_state * state );
int water_PH( double P, double H, water_state * state );
int water_PS( double P, double S, water_state * state );
int water_HS( double H, double S, water_state * state );
int water_TQ( double T, double Q, water_state * state );
int water_PQ( double P, double Q, water_state * state );

// Miscellaneous functions (will return -9.0e99 if the input is not valid).
double water_visc( double D, double T);
double water_cond( double D, double T);

namespace N_water_props
{
	const double T_crit = 647.096;
	const double P_crit = 22064.0;
	const double D_crit = 322.0;
	const double wmm = 18.015268;
	const double T_lower_limit = 273.2;
	const double T_upper_limit = 1500.0;
	const double P_lower_limit = 0.001;
	const double P_upper_limit = 50000.0;
	const double T_sat_min = 273.2;
	const double P_sat_min = 0.61343491;
	const double D_form_switch = 450.0;

	typedef struct
	{
		double x_low;
		double inv_dx;
		double y_low;
		double inv_dy;
		double c[6][6];
	}
	Element;

	void zero_state(water_state * __restrict state);

	void get_derivatives(
		const double x,
		const double y,
		const double D,
		const Element * __restrict e,
		double * __restrict f,
		double * __restrict dfdD,
		double * __restrict dfdD2,
		double * __restrict dfdT,
		double * __restrict dfdDdT,
		double * __restrict dfdT2
		);

	void get_two_phase_derivatives(
		const double x,
		const double y,
		const double D,
		const Element * __restrict e,
		double * __restrict f,
		double * __restrict dfdD,
		double * __restrict dfdT
		);

	void find_element(const double T, const double D, Element * __restrict e);

	void get_prop_derivatives(
		const double T,
		const double D,
		double * __restrict dPdD_T,
		double * __restrict dhdD_T,
		double * __restrict dsdD_T,
		double * __restrict dPdT_D,
		double * __restrict dhdT_D,
		double * __restrict dsdT_D,
		double * __restrict dDdP_T,
		double * __restrict dDdT_P,
		double * __restrict pres,
		double * __restrict enth,
		double * __restrict entr
		);
};

#endif

