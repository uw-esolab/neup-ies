#include "med_model.h"
#include "water_properties.h"
#include "saline.h"

//  #########################################################
int main()
{
    double 
        T = 315., 
        P = 0.007,
        S = 0.03;

    water_state ws; 
    water_TP(T, P, &ws);

    water_state bs;

    seawater_TPS(T, P, S, &bs);

    double h = bs.enth;

    return 0;
}

// ####################################################################
MED::MED(){};

MED::MED(int k)
{
    vapor_rate.clear();                          //Flow of the vapor rate for each n-effect
    feed_conc = .0335;                       //Average concentration of salt in seawater
    brine_rate.clear();                         //Brine flow rate 
    max_brine_conc = .067;                     //Maximum brine concentration
    this->k = k;                                   //Number of desired effects + 1 for starting values
    vbtemp.resize(k, 0.);            //Vapor and brine temperature vector creation
    vbtemp.at(k-1) = 40;                  //Fixing the final temperature of the brine and vapor
    enth_vapor.resize(k, 0.);            //Vapor enthalpy vector creation
    enth_brine.resize(k, 0.);            //Brine enthalpy vector creation
    enth_brine_nea.resize(k, 0.);       //allowance for flashing
    water_temp = 67.6;                     //Known starting/inlet temp of the DEMINERALIZED WATER [T2] (sCO2 outlet - delta_T_PCHE)
    //feed_temp  = 20                       //Known starting/inlet temp of the BRINE [T1]
    latentheat.resize(k, 0.);            //Latent heat of vapor vector creation
    tempchange = 3.;                          //Known temperature change, from Sharan paper (delta_T_NEA)
    water_rate = 308;                      //Is this the feed flow rate of the DEMINERALIZED WATER? m_dot_w2
    pressure.resize(k, 0.);                          //Pressure in MPa 
    brine_out.clear();
    nea = 3.; 

    /*
    NEA: upper value from MSF NEA is 1. The correlation Sharan used (JPME paper) is very difficult to interpret. 
    te difference is fairly small, perhaps 4% higher than Sharan with 4 effects. Setting NEA to 3 gives
    better agreement and could be due to looking at TD-TF in this paper, which is about 3C
    this gives almost perfect agreement with Sharan for 4 effects
    */

    std::vector<double> brine_conc_in;
    brine_conc_in.resize(k, feed_conc);
    
    double max_dif=1.;

    /*
    we had to assume brine concentrations to get the enthalpies. We iterate the brine concentrations (one loop does it)
    until we get convergence
    */
    
    
    while( max_dif > 0.01 )
    {
        //calculate Brine concentartion
        std::vector<double> brine_conc_out = iterate_MED(k,brine_conc_in); 
        
        //element-wise difference to previous iteration
        // max_dif = np.max(np.abs(np.divide(np.array(brine_conc_out),np.array(brine_conc_in))-1)) 
        for(size_t i=0; i<brine_conc_in.size(); i++)
        {
            double new_dif = brine_conc_out.at(i)/brine_conc_in.at(i)-1.;
            new_dif < 0 ? -new_dif : new_dif;

            if( new_dif > max_dif )
                max_dif = new_dif;
            
            brine_conc_in.at(i) = brine_conc_out.at(i);
        }
    }
}

// --------------------------------------------------------------------
std::vector<double> MED::iterate_MED(int k, std::vector<double>& brine_conc_in)
{

}

// --------------------------------------------------------------------
double MED::brine_flow_out(int i)
{
    double brine_out = feed_rate;
    for(size_t j=0; j<i+1; j++)
        brine_out += - vapor_rate.at(j); //Eq 6
    
    return brine_out;
}

// --------------------------------------------------------------------
double MED::brine_conctn(int i)
{
    double sum_vapor_rate = 0.;
    for(size_t j=0; j<i+1; j++)
        sum_vapor_rate += vapor_rate.at(j);

    double brine_concentration = (feed_rate*feed_conc)/(feed_rate-sum_vapor_rate);  //Eq 8
    return brine_concentration;
}

// --------------------------------------------------------------------
double MED::seaWaterSatH(double T, double S, double P)
{
    // the enthalpy at saturation pressure was pinging between gas and liquid. This clears it up
    water_state sw;
    seawater_TPS(T=T,P=P,S=S, &sw);
    double test_enth = sw.enth;
    if (sw.enth < 1000)
    {
        return sw.enth;
    }
    else
    {
        seawater_TPS(T-0.1, S, P, &sw);
        return sw.enth;
    }
}
