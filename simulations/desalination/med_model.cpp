#include "med_model.h"
#include "water_properties.h"
#include "saline.h"
#include "lp_lib.h"
#include <iostream>

//  #########################################################
int main()
{

    std::cout << "Enter the number of MED stages. A negative number will return all numbers in range 1...N:\n";
    int n_user = -1;
    std::cin >> n_user;

    if (n_user < 0)
    {
        for (int k = 1; k < -n_user; k++)
        {
            MED med_case(k);
            std::cout << k << "\t" << med_case.m3_per_day << "\n";
        }
    }
    else
    {
        MED med_case(n_user);
        std::cout << n_user << "\t" << med_case.m3_per_day << "\n";
    }

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
    m3_per_day = std::numeric_limits<double>::quiet_NaN();

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
        max_dif = 0.;

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

/*
Create a matrix handler structure
*/
struct _amat
{
protected:
    //internal data storage array
    double *_a;
    
public:
    int _nc;    //number of columns
    int _nr;    //number of rows
    
    _amat(){};  //empty initializer

    ~_amat()    //destructor
    {
        delete[] _a;    //free the data matrix on destruction
    };

    _amat(int nr, int nc)  //constructor
    {
        _nc = nc;
        _nr = nr;
        _a = new double[nr*nc];  //allocate memory
    };

    double& operator()(int i, int j)  //overload the () operator to mimic python's numpy convention
    {
        return _a[i*_nc+j];
    };

    double* array() //return address of internal data array
    {
        return _a;
    };
};

// --------------------------------------------------------------------
std::vector<double> MED::iterate_MED(int k, std::vector<double>& brine_conc_in)
{
    /* 
    Iterates to solve for MED states given a number of stages 'k' and brine concentrations. 

    Returns the new brine concentrations.
    */

    // step backwards through the stages to get the vapor temperatures
    for(int i=k-2; i>-1; i--)
        vbtemp.at(i) = vbtemp.at(i+1) + tempchange;
    
    // # A closer read of Sharan suggests this is the temperature of the first effect
    feed_temp = vbtemp.front();    
    
    water_state ws;
    double cp_water = std::numeric_limits<double>::quiet_NaN();

    for(int i=0; i<k; i++)
    {
        // evaluate enthalpy of saturated vapor
        water_TQ(vbtemp.at(i)+273.15, 1., &ws);
        enth_vapor.at(i) = ws.enth;

        pressure.at(i) = ws.pres/1000.;

        // brine enthalpy at liquid saturation
        enth_brine.at(i) = seaWaterSatH(vbtemp.at(i)+273.15, pressure.at(i), brine_conc_in.at(i));
        
        // brine enthalpy at condition subcooled by 'nea' degrees
        seawater_TPS(vbtemp.at(i) + 273.15 - nea, pressure.at(i), brine_conc_in.at(i), &ws);
        enth_brine_nea.at(i) = ws.enth;

        // evaluate enthalpy at saturated liquid
        water_TQ(vbtemp.at(i)+273.15, 0., &ws);

        // latent heat of evaporation for steam
        latentheat.at(i) = enth_vapor.at(i) -  ws.enth;

        //record the feed water specheat
        if (i == 0)
            cp_water = ws.cp;
    }
        
    double enth_feed = seaWaterSatH(feed_temp+273.15, pressure.back(), feed_conc);

    // Finding the vapor_rate for each n-effect. See ReTI/NE-2/Technical/Desalination Equations from Sharan for definitions of A and C matrices
    double* C = new double[k];
    // Create A matrix
    _amat A(k,k);

    A(0,0) = (-max_brine_conc/(max_brine_conc-feed_conc))*(enth_feed - enth_brine.at(0))+(enth_vapor.at(0)-enth_brine.at(0));
    
    C[0] = water_rate*(water_temp-(feed_temp+tempchange))* cp_water;

    //Creating the first row of the matrix
    for(int j=1; j<k; j++)           
    {    
        A(0,j) = (-max_brine_conc/(max_brine_conc-feed_conc))*(enth_feed - enth_brine[0]);
    }
                    
    for(int q=1; q<k; q++)            //Filling in the rest of the matrix
    {    
        //Below is the diagonal of the matrix
        A(q,q)= -(enth_vapor.at(q) - enth_brine.at(q)) + (max_brine_conc/(max_brine_conc-feed_conc))*(enth_brine_nea.at(q-1) - enth_brine.at(q));
        
        //Below is on column to the left of the diagonal of the matrix
        A(q,q-1)= latentheat.at(q-1)+((max_brine_conc/(max_brine_conc-feed_conc))-1)*(enth_brine_nea.at(q-1) - enth_brine.at(q));
        
        for(int j=0; j<q-1; j++) 
            //Below is everything to the left, under the q-1 diagonal. 
            A(q,j)= (max_brine_conc/(max_brine_conc-feed_conc)-1)*(enth_brine_nea.at(q-1) - enth_brine.at(q));
        for(int j=q+1; j<k; j++) 
            //Below is everything to the right of the q diagonal
            A(q,j)= (max_brine_conc/(max_brine_conc-feed_conc))*(enth_brine_nea.at(q-1) - enth_brine.at(q));
                    
        //For eqn 2-n, every term is related to a vapor rate
        C[q] = 0;
    }    

    // invA = np.linalg.inv(A);
    {
        double* vap_tmp = new double[k];
        // vapor_rate = np.dot(invA,C); //Solve for the vapor rate
        matinv(A.array(), C, k, vap_tmp);

        vapor_rate.resize(k);
        for (int i = 0; i < k; i++)
            vapor_rate.at(i) = vap_tmp[i];
        delete[] vap_tmp;
    }

    double distill = 0.;
    for (int i = 0; i < k; i++)
        distill += vapor_rate.at(i);
    
    //Calculate Brine feed flow with Eq 10
    feed_rate = distill*max_brine_conc/(max_brine_conc-feed_conc);

    std::vector<double> brine_conc; 
    for(int i=0; i<k; i++)
    {
        brine_rate.push_back(brine_flow_out(i));    //Updates brine_rate variable for every n-effect
        brine_conc.push_back(brine_conctn(i));
    }
    
    m3_per_day = distill*3600*24/1000.0; //from Sharan - convert the units from m3/s

    delete[] C; //free allocated memory

    return brine_conc;
}

// --------------------------------------------------------------------
double MED::brine_flow_out(int i)
{
    double brine_out = feed_rate;
    for(size_t j=0; j<(size_t)(i+1); j++)
        brine_out += - vapor_rate.at(j); //Eq 6
    
    return brine_out;
}

// --------------------------------------------------------------------
double MED::brine_conctn(int i)
{
    double sum_vapor_rate = 0.;
    for(size_t j=0; j< (size_t)(i + 1); j++)
        sum_vapor_rate += vapor_rate.at(j);
    
    return (feed_rate * feed_conc) / (feed_rate - sum_vapor_rate);  //Eq 8
}

// --------------------------------------------------------------------
double MED::seaWaterSatH(double T, double P, double S)
{
    // the enthalpy at saturation pressure was pinging between gas and liquid. This clears it up
    water_state sw;
    int ret = seawater_TPS(T, P, S, &sw, true);
    return sw.enth;
}

int matinv(double* A, double* b, int N, double* X)
{
    /*
    Perform matrix inversion to solve the linear system X=A\b

    Inputs:
        A   | double* - dimensions [N x N], row 1-> A[0], A[1], ... A[N]; row 2-> A[N+0], A[N+1], ... A[N+N-1], etc
        C   | double* - dimension [N], right-hand-side coefficients in linear equations
        N   | dimension of the square matrix

    Outputs:
        X   | double* - dimension [N], values solving the set of linear equations
    Returns:
        int - return code. 0 if successful. Other possible return codes are:
        --------------
        Solver status values
        UNKNOWNERROR            -5
        DATAIGNORED             -4
        NOBFP                   -3
        NOMEMORY                -2
        NOTRUN                  -1
        OPTIMAL                  0
        SUBOPTIMAL               1
        INFEASIBLE               2
        UNBOUNDED                3
        DEGENERATE               4
        NUMFAILURE               5
        USERABORT                6
        TIMEOUT                  7
        RUNNING                  8
        PRESOLVED                9


    This method uses LPSOLVE, which is an external open source linear algebra library. The source for this
    implementation was taken from NREL/SAM's SSC repository. https://github.com/nrel/ssc
    It is compiled as a .lib externally and linked during compilation of this program. The compiler 
    preprocessor definition should include "LPWINAPP" when compiling on windows.

    */
   //create the LP
    lprec* lp;

    /* Create a new LP model */
    lp = make_lp(0, N);
    if (lp == NULL) {
        fprintf(stderr, "Unable to create new LP model\n");
        return 0;
    }

    // automatically expand the matrix by adding rows
    set_add_rowmode(lp, TRUE);

    //storage for A matrix coefficients
    REAL* row = new REAL[(int)(N+1)];

    // set the objective function (all 1's - value doesn't matter), and set bounds to unconstrained
    for (int i = 0; i < N; i++)
    {
        row[i+1] = 1.;        //pathological objective function

        set_bounds(lp, i+1, -INFINITY, INFINITY);
    }
    set_obj_fn(lp, row);    //assign the objective function
    
    //add each row of the matrix
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            row[j+1] = A[i * N + j];         //row index is 1-based
        }
        add_constraint(lp, row, EQ, b[i]);   //set each row of the matrix, A[i]*x[i] = b[i]
    }

    set_maxim(lp);      //doesn't matter, choose to maximize

    set_add_rowmode(lp, FALSE);     //wrap up adding rows

    set_verbose(lp, CRITICAL);

    //write the problem to a file for debugging
    //char filename[100] = "C:/repositories/temp/lpsolve-dump.txt";
    //set_outputfile(lp, filename);
    //print_lp(lp);  //actually writes the file

    int ret = solve(lp);   //solves the problem
    //print_solution(lp, 1);  //prints the solution to file

    bool return_ok = ret == OPTIMAL || ret == SUBOPTIMAL;


    //retrieve variable values here
    if (return_ok)
    {
        get_variables(lp, X);
    }

    //clean up memory allocation
    delete[] row;

    //clean up memory allocation
    delete_lp(lp);
    return(0);
}