#include <math.h>

#include "../../nseos/crust.h"
#include "../../nseos/nuclear_en.h"
#include "../../nseos/coulomb.h"
#include "../../nseos/phyconst.h"
#include "../../nseos/mathconst.h"
#include "../../nseos/modeling.h"

double approximate_melting_temperature(struct compo comp, double nb_);
double eval_melting_temperature(
        struct parameters satdata, struct sf_params sparams, 
        double nb_);

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "ERROR: syntax is './ocrust_crystallization set.in'");
        return 1;
    }

    float b = 10.*log(2.);
    struct parameters satdata = assign_param(argv[1], b);
    double p = 3.0;
    struct sf_params sparams = fit_sf_params(satdata, p);

    double nb = 1.0e-8;

    while (nb < 2.5e-4)
    {
        fprintf(stdout, "%g %g\n", 
                nb, 
                eval_melting_temperature(satdata, sparams, nb)*1.16e10*1.e-9);

        nb += nb/10.;
    }

    return 0;
}

double approximate_melting_temperature(struct compo comp, double nb_)
{
    double gamma_m = 175.; // ion coupling parameter
    double zz = comp.aa*(1.-comp.del)/2.;

    return  zz*zz*ALPHAFS*HBARC/gamma_m*pow(4.*PI*nb_/3./comp.aa, 1./3.);
}

double eval_melting_temperature(
        struct parameters satdata, struct sf_params sparams,
        double nb_)
{
        double guess_liq[3] = {60., 0.10, 0.1595};
        struct compo comp_liq = calc_ocrust3d_composition(nb_, 0., 
                "sol", guess_liq, satdata, sparams);

        double tt_init = approximate_melting_temperature(comp_liq, nb_)*2.;

        double fws_sol, fws_liq;

        double tt = tt_init; 

        do
        {
            comp_liq = calc_ocrust3d_composition(nb_, tt, 
                    "liq", guess_liq, satdata, sparams); 

            fws_sol = calc_crust_ws_cell_free_energy_density(
                    satdata, sparams, comp_liq, nb_, tt, "sol");
            fws_liq = calc_crust_ws_cell_free_energy_density(
                    satdata, sparams, comp_liq, nb_, tt, "liq");

            tt -= 0.0001;

            if (tt <= 0.)
            {
                fprintf(stderr, "ERROR: sign of T_m cannot be negative!\n");
                return -1;
            }
        }
        while(fws_liq < fws_sol);

        return tt;
}
