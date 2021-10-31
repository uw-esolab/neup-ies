#include "saline.h"
#include <exception>
#include <math.h>
#include <limits>

int seawater_TPS(double T, double P, double S, water_state * brine_state)
{
    /*
    Model seawater with standard IAPWS-08

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [MPa]
    S : float
        Salinity, [kg/kg]

    Returns
    -------
    rho : float
        Density, [kg/m³]
    h : float
        Specific enthalpy, [kJ/kg]
    s : float
        Specific entropy, [kJ/kg·K]
    cp : float
        Specific isobaric heat capacity, [kJ/kg·K]
    
    */

    // Constants
    double _g, _gt, _gp, _gtt, _gtp, _gpp, _gs, _gsp;
    
    _g=0.;
    _gt=0.;
    _gp=0.;
    _gtt=0.;
    _gtp=0.;
    _gpp=0.;
    _gs=0.;
    _gsp=0.;

    // ----------------------------
    double Rm;
    double Sn;
    double S_;
    double Ms;
    double Po;
    double To;
    Rm = 8.314472;
    Sn = 0.03516504;
    S_ = Sn*40/35;
    Ms = 31.4038218;
    Po = 0.101325;
    To = 273.15;

    // Tables and data
    int I[64] = 
        {1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 2, 3, 4, 2, 3, 4, 2, 3, 4,
        2, 4, 2, 2, 3, 4, 5, 2, 3, 4, 2, 3, 2, 3, 2, 3, 2, 3, 4, 2, 3, 2,
        3, 2, 2, 2, 3, 4, 2, 3, 2, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2};
    int J[64] = 
        {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
        5, 5, 6, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 0, 1, 1, 2,
        2, 3, 4, 0, 0, 0, 1, 1, 2, 2, 3, 4, 0, 0, 1, 2, 3, 0, 1, 2};
    int K[64] = 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5};
    double G[64] =
    {
        0.581281456626732e4, 0.141627648484197e4, -0.243214662381794e4,
        0.202580115603697e4, -0.109166841042967e4, 0.374601237877840e3,
        -0.485891069025409e2, 0.851226734946706e3, 0.168072408311545e3,
        -0.493407510141682e3, 0.543835333000098e3, -0.196028306689776e3,
        0.367571622995805e2, 0.880031352997204e3, -0.430664675978042e2,
        -0.685572509204491e2, -0.225267649263401e3, -0.100227370861875e2,
        0.493667694856254e2, 0.914260447751259e2, 0.875600661808945,
        -0.171397577419788e2, -0.216603240875311e2, 0.249697009569508e1,
        0.213016970847183e1, -0.331049154044839e4, 0.199459603073901e3,
        -0.547919133532887e2, 0.360284195611086e2, 0.729116529735046e3,
        -0.175292041186547e3, -0.226683558512829e2, -0.860764303783977e3,
        0.383058066002476e3, 0.694244814133268e3, -0.460319931801257e3,
        -0.297728741987187e3, 0.234565187611355e3, 0.384794152978599e3,
        -0.522940909281335e2, -0.408193978912261e1, -0.343956902961561e3,
        0.831923927801819e2, 0.337409530269367e3, -0.541917262517112e2,
        -0.204889641964903e3, 0.747261411387560e2, -0.965324320107458e2,
        0.680444942726459e2, -0.301755111971161e2, 0.124687671116248e3,
        -0.294830643494290e2, -0.178314556207638e3, 0.256398487389914e2,
        0.113561697840594e3, -0.364872919001588e2, 0.158408172766824e2,
        -0.341251932441282e1, -0.316569643860730e2, 0.442040358308000e2,
        -0.111282734326413e2, -0.262480156590992e1, 0.704658803315449e1,
        -0.792001547211682e1
    };

    // # Constants
    // Check input in range of validity
    if ((T <= 261) || (T > 353) || (P <= 0) || (P > 100) || (S < 0) || (S > 0.12))
        //  warnings.warn("Incoming out of bound")s
        // throw std::runtime_error("Incoming parameter out of bounds");
        return -1;

    // Get the state of pure water
    water_state w_state;
    water_TP(T, P*1000., &w_state);
    water_state w_ref;
    water_TP(To, Po*1000., &w_ref);

    // Calculate the state of the salt component
    S_ = Sn*40/35;
    double X = pow((S/S_),0.5);
    double tau = (T-To)/40;
    double pi = (P-Po)/100;

    double m = S/(1-S)/Ms;

    _g=0.;
    _gt=0.;
    _gp=0.;
    _gtt=0.;
    _gtp=0.;
    _gpp=0.;
    _gs=0.;
    _gsp=0.;

    // # Calculate only for some salinity
    if (S != 0)
    {
        for(int n=0; n<64; n++)
        {
            double i = I[n];
            double j = J[n];
            double k = K[n];
            double gi = G[n];

            if (i == 1)
            {
                _g += gi*X*X*log(X)*pow(tau,j)*pow(pi,k);
                _gs += gi*(2*log(X)+1)*pow(tau,j)*pow(pi,k);
            }
            else
            {
                _g += gi*pow(X,i)*pow(tau,j)*pow(pi,k);
                _gs += i*gi*pow(X,(i-2))*pow(tau,j)*pow(pi,k);
            }
            if (j >= 1)
            {
                if (i == 1)
                    _gt += gi*pow(X,2)*log(X)*j*pow(tau,(j-1))*pow(pi,k);
                else
                    _gt += gi*pow(X,i)*j*pow(tau,(j-1))*pow(pi,k);
            }
            if (k >= 1)
            {
                _gp += k*gi*pow(X,i)*pow(tau,j)*pow(pi,(k-1));
                _gsp += i*k*gi*pow(X,(i-2))*pow(tau,j)*pow(pi,(k-1));
            }
            if (j >= 2)
                _gtt += j*(j-1)*gi*pow(X,i)*pow(tau,(j-2))*pow(pi,k);
            if ((j >= 1) && (k >= 1))
                _gtp += j*k*gi*pow(X,i)*pow(tau,(j-1))*pow(pi,(k-1));
            if (k >= 2)
                _gpp += k*(k-1)*gi*pow(X,i)*pow(tau,j)*pow(pi,(k-2));
        }
    }

    // Conversions
    _g = _g*1e-3;               //gibbs -> h - T*s
    _gt = _gt/40*1e-3;          //d(gibbs)/dT
    _gp = _gp/100*1e-6;         //d(gibbs)/dp = spec vol.
    _gtt = _gtt/1600.*1e-3;
    _gtp = _gtp/40/100*1e-6;
    _gpp = _gpp/10000*1e-6;
    _gs = _gs/S_/2*1e-3;
    _gsp = _gsp/S_/2/100*1e-6;


    // Calculate combined values
    brine_state->temp = T;
    brine_state->pres = P;
    
    
    // Gibbs free energy of pure water
    double gfe_water = w_state.enth - T * w_state.entr;
    // Gibbs combined
    double gfe_brine = gfe_water + _g;
    // entropy combined
    double entr_brine = w_state.entr - _gt;

    //enthalpy combined
    // self.s = -prop["gt"]
    brine_state->entr = entr_brine;
    // self.h = prop["g"]-T*prop["gt"]
    brine_state->enth = gfe_brine + T * entr_brine;


    // self.v = prop["gp"]
    // self.rho = 1./prop["gp"]
    brine_state->dens = 1./(_gp + 1./w_state.dens);

    // self.cp = -T*prop["gtt"]
    brine_state->cp = -T*(_gtt - w_state.cp/T);
    // self.cv = T*(prop["gtp"]**2/prop["gpp"]-prop["gtt"])
    // self.u = prop["g"]-T*prop["gt"]-P*1000*prop["gp"]
    // self.a = prop["g"]-P*1000*prop["gp"]
    // self.alfav = prop["gtp"]/prop["gp"]
    // self.betas = -prop["gtp"]/prop["gtt"]
    // self.xkappa = -prop["gpp"]/prop["gp"]
    // self.ks = (prop["gtp"]**2-prop["gtt"]*prop["gpp"])/prop["gp"] / \
    //     prop["gtt"]
    // self.w = prop["gp"]*(prop["gtt"]*1000/(prop["gtp"]**2 -
    //                      prop["gtt"]*1000*prop["gpp"]*1e-6))**0.5

    return 1;
}
