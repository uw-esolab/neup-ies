#include "water_properties.h"
#include "saline.h"

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