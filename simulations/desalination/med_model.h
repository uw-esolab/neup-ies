/*
MED Model 
ReTI Research Lab / ESOLab 
Created by Grace Stanke, Ben Lindley, Mike Wagner

Based on Sharan 2019
As with Sharan, gives optimal water output with 4 effects
*/

#include <vector>

int matinv(double* A, double* b, int n, double* X);

double water_value_in_arizona()
{
    // https://efc.web.unc.edu/2014/09/23/conservation-water-rates-arizona-utilities-using-rates-discourage-discretionary-water-use/
    double price_per_1000_gallons=3.5;
    double gallons_per_m3 = 264.172;
    return price_per_1000_gallons/1000*gallons_per_m3;
}

class MED
{
    /*
    values we ultimately expect to pass to the MED model:
    k=4
    water rate - solved through EES mass balance
    water_temp - solved through EES mass and temperature balance
    */

    int k;                                   //Number of desired effects + 1 for starting values
    double feed_conc;                       //Average concentration of salt in seawater
    double max_brine_conc;                     //Maximum brine concentration
    double water_temp;                     //Known starting/inlet temp of the DEMINERALIZED WATER [T2] (sCO2 outlet - delta_T_PCHE)
    double tempchange;                          //Known temperature change, from Sharan paper (delta_T_NEA)
    double water_rate;                      //Is this the feed flow rate of the DEMINERALIZED WATER? m_dot_w2
    double feed_rate;
    double feed_temp;
    double nea; 
    std::vector<double> vapor_rate;                          //Flow of the vapor rate for each n-effect
    std::vector<double> brine_rate;                         //Brine flow rate 
    std::vector<double> vbtemp;            //Vapor and brine temperature vector creation
    std::vector<double> enth_vapor;            //Vapor enthalpy vector creation
    std::vector<double> enth_brine;            //Brine enthalpy vector creation
    std::vector<double> enth_brine_nea;       //allowance for flashing
    std::vector<double> latentheat;            //Latent heat of vapor vector creation
    std::vector<double> pressure;                          //Pressure in MPa 
    std::vector<double> brine_out;

public:
    MED();
    MED(int k);
    std::vector<double> iterate_MED(int k, std::vector<double> &brine_conc_in);
    double brine_flow_out(int i);
    double brine_conctn(int i);
    double seaWaterSatH(double T, double P, double S);
    double m3_per_day;
};