#include <math.h>

#include "../../nseos/crust.h"
#include "../../nseos/nuclear_en.h"
#include "../../nseos/coulomb.h"
#include "../../nseos/phyconst.h"
#include "../../nseos/mathconst.h"
#include "../../nseos/modeling.h"

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        fprintf(stderr, 
                "ERROR: syntax is './ocrust_crystallization "
                "set.in ocrust.out eos.out'\n");
        return 1;
    }

    // nuclear modeling
    float b = 10.*log(2.);
    struct parameters satdata = assign_param(argv[1], b);
    double p = 3.0;
    struct sf_params sparams = fit_sf_params(satdata, p);

    double nb = 1.0e-8; // initial baryon density
    double tm;
    struct compo comp;

    FILE *ocrust = fopen(argv[2],"w+");
    FILE *eos = fopen(argv[3],"w+");

    while (nb < 2.5e-4) // up to neutron drip density
    {
        tm = eval_melting_temperature(satdata, sparams, nb, &comp);

        print_state_crust(satdata, sparams, comp, nb, tm, "sol", ocrust, eos);

        nb += nb/10.;
    }

    fclose(ocrust);

    return 0;
}
